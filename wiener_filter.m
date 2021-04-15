function wiener_filter(s,x1,x2,x3)
    R=xcorr(x1,x1);
    %rxd=xcorr(x1,conj(s));
    %w=inv(R)*rxd

end