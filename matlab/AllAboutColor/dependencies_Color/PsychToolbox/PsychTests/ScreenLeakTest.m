function screenleaktest% ScreenLeakTest measures memory leakage associated with opening and % closing windows.clear allbytes;fprintf('%8ld free, %8ld temp free. After initial CLEAR ALL.\n',bytes,bytes('TempFree'));screen('windows');fprintf('%8ld free, %8ld temp free. After SCREEN(''windows'').\n',bytes,bytes('TempFree'));% whichScreen=0;%window=QT('OpenWindow',whichScreen);SCREEN(0,'OpenWindow');SCREEN(0,'OpenOffscreenWindow');fprintf('%8ld free, %8ld temp free. After SCREEN(''OpenWindow'').\n',bytes,bytes('TempFree'));%QT('CloseWindow',window);SCREEN('CloseAll');fprintf('%8ld free, %8ld temp free. After SCREEN(''CloseAll'').\n',bytes,bytes('TempFree'));%screen close?clear allbytes;fprintf('%8ld free, %8ld temp free. After CLEAR ALL.\n',bytes,bytes('TempFree'));