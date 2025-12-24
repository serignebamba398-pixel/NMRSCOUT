from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from PIL import Image


@dataclass(frozen=True)
class Options:
    sizes: tuple[int, ...] = (16, 24, 32, 48, 64, 128, 256)
    corner_threshold: int = 250  # consider "near white" background
    corner_tolerance: int = 40   # allow JPEG-ish off-white in background
    edge_clear_px: int = 12      # also clear near-white pixels near borders
    pad_ratio: float = 0.0       # keep no empty border around logo


def _is_near_white(rgb: tuple[int, int, int], threshold: int, tolerance: int) -> bool:
    r, g, b = rgb
    return r >= threshold - tolerance and g >= threshold - tolerance and b >= threshold - tolerance


def _flood_fill_remove_corner_bg(img: Image.Image, threshold: int, tolerance: int) -> Image.Image:
    """Turn corner-connected near-white background pixels transparent.

    Works well for JPEG sources with a white matte background.
    """

    rgba = img.convert("RGBA")
    px = rgba.load()
    w, h = rgba.size

    visited = [[False] * w for _ in range(h)]

    def enqueue(queue: list[tuple[int, int]], x: int, y: int) -> None:
        if x < 0 or y < 0 or x >= w or y >= h:
            return
        if visited[y][x]:
            return
        visited[y][x] = True
        queue.append((x, y))

    queue: list[tuple[int, int]] = []
    # Seed from corners
    for x, y in [(0, 0), (w - 1, 0), (0, h - 1), (w - 1, h - 1)]:
        enqueue(queue, x, y)

    while queue:
        x, y = queue.pop()
        r, g, b, a = px[x, y]
        if a == 0:
            # already transparent, continue flood
            for nx, ny in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                enqueue(queue, nx, ny)
            continue

        if _is_near_white((r, g, b), threshold, tolerance):
            px[x, y] = (r, g, b, 0)
            for nx, ny in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                enqueue(queue, nx, ny)

    return rgba


def _clear_near_white_edges(rgba: Image.Image, threshold: int, tolerance: int, edge_px: int) -> Image.Image:
    """Clear near-white pixels in an outer band around the image.

    Useful for JPEG sources where background isn't perfectly connected in flood-fill
    due to compression artifacts.
    """

    if edge_px <= 0:
        return rgba

    img = rgba.copy().convert("RGBA")
    px = img.load()
    w, h = img.size
    edge_px = min(edge_px, (min(w, h) // 2) if min(w, h) > 1 else 1)

    def in_edge_band(x: int, y: int) -> bool:
        return x < edge_px or y < edge_px or x >= (w - edge_px) or y >= (h - edge_px)

    for y in range(h):
        for x in range(w):
            if not in_edge_band(x, y):
                continue
            r, g, b, a = px[x, y]
            if a == 0:
                continue
            if _is_near_white((r, g, b), threshold, tolerance):
                px[x, y] = (r, g, b, 0)

    return img


def _dehalo_white_matte(rgba: Image.Image) -> Image.Image:
    """Remove white fringe from anti-aliased edges (assumes white matte)."""

    img = rgba.copy().convert("RGBA")
    px = img.load()
    w, h = img.size

    for y in range(h):
        for x in range(w):
            r, g, b, a = px[x, y]
            if a == 0 or a == 255:
                continue

            # Solve fg from: observed = fg*a + white*(1-a)
            af = a / 255.0
            inv = (1.0 - af)
            # observed and matte are in [0,255]
            # fg = (observed - matte*(1-a)) / a
            def solve(ch: int) -> int:
                fg = (ch - 255.0 * inv) / af
                return max(0, min(255, int(round(fg))))

            nr = solve(r)
            ng = solve(g)
            nb = solve(b)
            px[x, y] = (nr, ng, nb, a)

    return img


def _tight_bbox_alpha(img: Image.Image) -> tuple[int, int, int, int] | None:
    """Return bbox of non-transparent pixels, or None if empty."""

    if img.mode != "RGBA":
        img = img.convert("RGBA")
    alpha = img.split()[-1]
    return alpha.getbbox()


def _scale_to_fill(img: Image.Image, pad_ratio: float) -> Image.Image:
    """Crop to content and scale to *cover* a square canvas.

    This avoids empty left/right borders for non-square logos by scaling up
    and center-cropping (like CSS background-size: cover).
    """

    rgba = img.convert("RGBA")
    bbox = _tight_bbox_alpha(rgba)
    if not bbox:
        return rgba

    content = rgba.crop(bbox)
    cw, ch = content.size
    side = max(cw, ch)

    pad_ratio = max(0.0, min(0.25, float(pad_ratio)))
    inner = int(round(side * (1.0 - 2.0 * pad_ratio)))
    inner = max(1, inner)

    # Scale content to cover inner x inner
    scale = max(inner / cw, inner / ch)
    nw = max(1, int(round(cw * scale)))
    nh = max(1, int(round(ch * scale)))
    scaled = content.resize((nw, nh), resample=Image.Resampling.LANCZOS)

    # Center-crop to inner square
    left = max(0, (nw - inner) // 2)
    top = max(0, (nh - inner) // 2)
    cropped = scaled.crop((left, top, left + inner, top + inner))

    out = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    px = (side - inner) // 2
    py = (side - inner) // 2
    out.paste(cropped, (px, py))
    return out


def _ensure_square(img: Image.Image) -> Image.Image:
    rgba = img.convert("RGBA")
    w, h = rgba.size
    if w == h:
        return rgba
    side = max(w, h)
    out = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    out.paste(rgba, ((side - w) // 2, (side - h) // 2))
    return out


def generate_ico(source: Path, target: Path, options: Options) -> None:
    base = Image.open(source)

    # Corner-based background removal + dehalo, so JPEG-with-white can become clean RGBA
    rgba = _flood_fill_remove_corner_bg(
        base, threshold=options.corner_threshold, tolerance=options.corner_tolerance
    )
    rgba = _clear_near_white_edges(
        rgba,
        threshold=options.corner_threshold,
        tolerance=options.corner_tolerance,
        edge_px=options.edge_clear_px,
    )
    rgba = _dehalo_white_matte(rgba)
    rgba = _scale_to_fill(rgba, pad_ratio=options.pad_ratio)
    rgba = _ensure_square(rgba)

    sizes = [s for s in options.sizes if s > 0]
    # Pillow writes a valid ICO when given multiple sizes.
    tmp = target.with_suffix(target.suffix + ".tmp")
    tmp.parent.mkdir(parents=True, exist_ok=True)

    rgba.save(tmp, format="ICO", sizes=[(s, s) for s in sizes])
    os.replace(tmp, target)


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    source = repo / "assets" / "nmrscout.png"
    target = repo / "assets" / "icon.ico"

    if not source.exists():
        print(f"ERROR: source not found: {source}")
        return 2

    options = Options()
    generate_ico(source, target, options)
    print(f"Wrote {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
