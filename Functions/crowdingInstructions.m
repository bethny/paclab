function crowdingInstructions(currBlock, w)

%==============================================================
% Adapted from Dan's displayInstructions.m
%==============================================================
grey = GrayIndex(w,0.5);
Screen('FillRect', w, grey);
Screen('TextSize', w, []);

instText = [];

%Practice Instructions
if currBlock == 0 
    instText = ['Press the LEFT ARROW KEY if the center bar\n'...
        'is tilted to the left (counterclockwise).\n'...
        'Press the RIGHT ARROW KEY if the center bar\n'...
        'is tilted to the right (clockwise).\n'...
        'If you respond correctly, you will hear a beep. \n\n'...
        'If you respond incorrectly, you will hear a lower beep. \n\n\n'...  
        'If you need to exit the trial, press q. \n\n'...
        'We will start with some practice trials. \n\n'...
        'Press the RIGHT ARROW KEY to continue.'];
elseif currBlock == 1
    instText = ['The practice block is now complete.\n\n'...
                'There will be a number of experimental blocks with breaks in between.\n\n'...                
                'Press the RIGHT ARROW KEY to continue.\n\n'];
else
    instText = 'Press the RIGHT ARROW KEY to continue.\n\n';
end
DrawFormattedText(w, instText, 'Center', 'Center', [255 255 255]);

Screen('Flip', w);
FlushEvents('keyDown');
GetChar;
Screen('TextSize', w, 20);
Screen('Flip', w);

