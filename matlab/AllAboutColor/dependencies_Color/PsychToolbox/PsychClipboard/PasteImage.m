function M = PasteImage% M = PasteImage% % Paste the clipboard into matrix M.  The clipboard should be a PICT. % M comes out as an image.  If the PICT had three color channels, this% information gets averaged, so that the returned image is grayscale.% Image data should be scaled between 0 and 1. After multiplying by % 255 it should be suitable for SCREEN 'PutImage'.% % BUG & WORK AROUND:% % PasteImage fails with some, but not all, images copied to the% clipboard from Adobe Photoshop. Copying the image from the clipboard% into another application (e.g. Microsoft Word) and then from that% application to the clipboard, seems to fix the problem. (Guessing as% to the cause, it's possible that Photoshop is embedding the% picture's ICC profile into the picture's comment, making the picture% odd in a way that PasteImage ought to but doesn't cope with, and% that Word is stripping the profile out. Embedding the ICC profile is% generally considered a good thing, enhancing the reproducibility of% the image, but, if this diagnosis is right, PasteImage isn't coping% with it.) Reported by Geoff Loftus 8/1/00.% % Also see PasteCImage, CopyImage, SCREEN.% 8/19/94		dhb		Link in new pixmap handling code.% 8/4/00		dgp		Mention Photoshop-related bug, and tentative diagnosis.