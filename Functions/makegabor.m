function m=makegabor(gaborDimPix,freq,sigma,ori,VgaborDimPix,transformOut)
    
    res=[round(gaborDimPix) round(gaborDimPix) ];
    phase = 180*round(rand*2);%*CoinFlip(1,.5)-1;
    sc = sigma;
    tilt = 90+ori;
    x=res(1)/2;
    y=res(2)/2;
    r=VgaborDimPix/2;
    sf = freq;

    [gab_x gab_y] = meshgrid(-res(1)/2:(res(1)/2-1), -res(2)/2:(res(2)/2-1));
    a=cos((pi*tilt/180))*sf*360;
    b=sin((pi*tilt/180))*sf*360;

    x_factor=-1*(gab_x).^2;
    y_factor=-1*(gab_y).^2;
    sinWave = sin(pi/180*(a*(gab_x - x) + b*(gab_y - y)+phase));

    for xx=-res(1)/2:(res(1)/2-1)
     for yy=-res(2)/2:(res(2)/2-1)
         if xx^2+yy^2 > r^2
             sinWave(xx+res(1)/2+1,yy+res(2)/2+1)=0;
         end
     end
    end

    varScale=2*sc^2;
    m_tmp = (exp(x_factor/varScale+y_factor/varScale).*sinWave)';
    if transformOut
        m = m_tmp;
    else
        m = m_tmp*0.5*127.5 + 127.5;
    end
    
    save('makegabor.mat')