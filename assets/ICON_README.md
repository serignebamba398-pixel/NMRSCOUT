# App icon (Windows)

This project is packaged with `electron-builder`.

## What to add

- Create a Windows icon file: `assets/icon.ico`
- Recommended sizes inside the `.ico`: 16, 24, 32, 48, 64, 128, 256 px

Tip: you typically author the logo at a larger size (e.g. 512x512 or 1024x1024) and then export an `.ico` that embeds the sizes above. Windows will pick the closest size depending on where the icon is shown (taskbar, Start menu, Explorer, shortcuts).

## How to enable it

Once `assets/icon.ico` exists, ensure the build config in the root `package.json` contains:

```json
{
  "build": {
    "icon": "assets/icon.ico"
  }
}
```

Then rebuild:

- `npm run pre-dist`
- `npm run dist`

Tip: if you only have a PNG logo, export a square transparent PNG (>= 512x512) and we can generate `assets/icon.ico` from it.
