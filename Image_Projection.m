        vT = 10;
        img_orig = importdata('lena_full.png');
        img = img_orig/(vT+1);
        for i = 1:vT
            img = img + imtranslate(img_orig/(vT+1),[i,0]);
        end
        figure;
        imshow(img_orig);
        figure;
        imshow(img);
        imwrite(img,'lena_blurry.png');