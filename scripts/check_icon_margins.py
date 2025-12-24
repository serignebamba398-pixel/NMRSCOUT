from __future__ import annotations

from pathlib import Path

from PIL import Image


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    ico = repo / "assets" / "icon.ico"
    if not ico.exists():
        print(f"Missing: {ico}")
        return 2

    im = Image.open(ico)
    frames = []
    try:
        while True:
            frames.append(im.copy().convert("RGBA"))
            im.seek(im.tell() + 1)
    except EOFError:
        pass

    if not frames:
        print("No frames found in ICO")
        return 3

    largest = max(frames, key=lambda f: f.size[0])
    w, h = largest.size
    alpha = largest.split()[-1]
    bbox = alpha.getbbox()

    print(f"frames={len(frames)}")
    print(f"largest={w}x{h}")
    print(f"alpha_bbox={bbox}")

    if not bbox:
        print("All transparent")
        return 4

    l, t, r, b = bbox
    print(f"margins: left={l} top={t} right={w - r} bottom={h - b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
