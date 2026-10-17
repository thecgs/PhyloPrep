import { readFile } from 'node:fs/promises';
import vm from 'node:vm';
const context = vm.createContext({ Blob });
vm.runInContext(await readFile('site/browser-save.js', 'utf8'), context);
export const PhyloPrepSave = context.PhyloPrepSave;
