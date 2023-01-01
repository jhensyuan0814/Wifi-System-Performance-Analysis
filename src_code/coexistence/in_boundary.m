function new_MS = in_boundary(MS, side)
    new_MS = MS;
    for i = 1: length(MS)
        if abs(real(MS)) > side / 2
            new_MS(i) = (side - abs(real(MS))) * sign(real(MS)) + 1i * imag(MS);
        end
        if abs(imag(MS)) > side / 2
            new_MS(i) = real(MS) + 1i * (side - abs(imag(MS))) * sign(imag(MS));
        end
    end
end