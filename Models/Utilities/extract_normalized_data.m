function [s1, o1, s2, o2] = extract_normalized_data(data, t)
    s1 = data(t, 3) / 10;
    o1 = data(t, 4) / 10;
    s2 = data(t, 5) / 10;
    o2 = data(t, 6) / 10;
end