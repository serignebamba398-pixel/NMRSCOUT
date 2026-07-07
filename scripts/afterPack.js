/**
 * afterPack hook for electron-builder
 * Cleans up Windows reserved filename "nul" that gets created during asar unpacking
 */
const fs = require('fs');
const path = require('path');

exports.default = async function(context) {
    const appOutDir = context.appOutDir;

    // Path where the nul file might be created
    const nulPath = path.join(appOutDir, 'resources', 'app.asar.unpacked', 'backend', 'nul');

    // Try to delete using rimraf-style approach for Windows reserved names
    try {
        if (fs.existsSync(nulPath)) {
            console.log('[afterPack] Found nul file, attempting to delete:', nulPath);

            // Use the \\?\ prefix for Windows reserved names
            const uncPath = '\\\\?\\' + nulPath;
            fs.unlinkSync(uncPath);
            console.log('[afterPack] Successfully deleted nul file');
        }
    } catch (err) {
        // If direct deletion fails, try alternative methods
        console.warn('[afterPack] Could not delete nul file directly:', err.message);

        try {
            // Try using spawn to call del command
            const { execSync } = require('child_process');
            execSync(`del /f /q "\\\\?\\${nulPath}"`, { stdio: 'inherit' });
            console.log('[afterPack] Deleted nul file via cmd');
        } catch (cmdErr) {
            console.warn('[afterPack] cmd deletion also failed:', cmdErr.message);
        }
    }
};
