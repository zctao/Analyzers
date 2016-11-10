#./parser input output
sed -E '/\*{2,}/d;
        /\*\ +(Row)/d;
        s/\*\ +//;
        s/\ \*\ */:/g;
        s/^[0-9]+://;
        s/:$//'< $1 > $2