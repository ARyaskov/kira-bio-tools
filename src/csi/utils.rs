pub fn reg2bin(start: usize, end: usize, min_shift: u8, depth: u8) -> usize {
    let end = end.saturating_sub(1);
    let mut s = min_shift as usize;
    let mut t = ((1 << (depth as usize * 3)) - 1) / 7;

    for _ in 0..depth {
        if start >> s == end >> s {
            return t + (start >> s);
        }
        s += 3;
        t = (t - 1) / 8;
    }

    0
}
