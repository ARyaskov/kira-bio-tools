use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AnnotateMode {
    pub replace_all: bool,
    pub replace_missing: bool,
    pub replace_non_missing: bool,
    pub set_or_append: bool,
    pub carry_over_missing: bool,
    pub match_value: bool,
}

impl AnnotateMode {
    pub fn default_mode() -> Self {
        Self {
            replace_all: true,
            replace_missing: false,
            replace_non_missing: false,
            set_or_append: false,
            carry_over_missing: false,
            match_value: false,
        }
    }

    pub fn parse(spec: &str) -> (Self, &str) {
        let mut mode = Self::default();
        let mut idx = 0;
        let chars: Vec<char> = spec.chars().collect();

        while idx < chars.len() {
            match chars[idx] {
                '+' => {
                    mode.replace_all = false;
                    mode.replace_missing = true;
                }
                '-' => {
                    mode.replace_all = false;
                    mode.replace_non_missing = true;
                }
                '=' => {
                    mode.replace_all = false;
                    mode.set_or_append = true;
                }
                '.' => {
                    mode.carry_over_missing = true;
                }
                '~' => {
                    mode.match_value = true;
                }
                _ => break,
            }
            idx += 1;
        }

        if !mode.replace_missing
            && !mode.replace_non_missing
            && !mode.set_or_append
            && !mode.match_value
        {
            mode.replace_all = true;
        }

        (mode, &spec[idx..])
    }

    pub fn should_transfer(
        &self,
        src_is_missing: bool,
        dst_exists: bool,
        dst_is_missing: bool,
    ) -> bool {
        if self.match_value {
            return false;
        }

        if self.set_or_append {
            return !src_is_missing || self.carry_over_missing;
        }

        if self.replace_non_missing {
            return dst_exists && !dst_is_missing && (!src_is_missing || self.carry_over_missing);
        }

        if self.replace_missing {
            return (!dst_exists || dst_is_missing) && (!src_is_missing || self.carry_over_missing);
        }

        if self.replace_all {
            return !src_is_missing || self.carry_over_missing;
        }

        false
    }

    pub fn should_append(&self) -> bool {
        self.set_or_append
    }
}

impl fmt::Display for AnnotateMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut prefix = String::new();
        if self.carry_over_missing {
            prefix.push('.');
        }
        if self.replace_missing {
            prefix.push('+');
        }
        if self.replace_non_missing {
            prefix.push('-');
        }
        if self.set_or_append {
            prefix.push('=');
        }
        if self.match_value {
            prefix.push('~');
        }
        write!(f, "{}", prefix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_modes() {
        let (mode, rest) = AnnotateMode::parse("TAG");
        assert!(mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("+TAG");
        assert!(mode.replace_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("-TAG");
        assert!(mode.replace_non_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("=TAG");
        assert!(mode.set_or_append);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse(".TAG");
        assert!(mode.carry_over_missing);
        assert!(mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse(".+TAG");
        assert!(mode.carry_over_missing);
        assert!(mode.replace_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("~ID");
        assert!(mode.match_value);
        assert_eq!(rest, "ID");
    }

    #[test]
    fn test_should_transfer() {
        let mode_tag = AnnotateMode::default_mode();
        assert!(mode_tag.should_transfer(false, false, false));
        assert!(mode_tag.should_transfer(false, true, false));
        assert!(mode_tag.should_transfer(false, true, true));
        assert!(!mode_tag.should_transfer(true, false, false));

        let (mode_plus, _) = AnnotateMode::parse("+TAG");
        assert!(mode_plus.should_transfer(false, false, false));
        assert!(!mode_plus.should_transfer(false, true, false));
        assert!(mode_plus.should_transfer(false, true, true));

        let (mode_minus, _) = AnnotateMode::parse("-TAG");
        assert!(!mode_minus.should_transfer(false, false, false));
        assert!(mode_minus.should_transfer(false, true, false));
        assert!(!mode_minus.should_transfer(false, true, true));

        let (mode_dot, _) = AnnotateMode::parse(".TAG");
        assert!(mode_dot.should_transfer(true, false, false));
        assert!(mode_dot.should_transfer(true, true, false));
    }
}
