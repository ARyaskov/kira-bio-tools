use super::super::structs::annotate_mode::AnnotateMode;
use super::field_metadata::is_missing_value;

pub fn merge_append_values(existing: &str, new_val: &str) -> String {
    if existing.is_empty() || existing == "." {
        new_val.to_string()
    } else if new_val.is_empty() || new_val == "." {
        existing.to_string()
    } else {
        format!("{},{}", existing, new_val)
    }
}

pub fn merge_default_values(
    mode: AnnotateMode,
    vcf_value: Option<&str>,
    db_value: Option<&str>,
) -> String {
    let vcf_val = vcf_value.unwrap_or(".");
    let db_val = db_value.unwrap_or(".");

    if mode.set_or_append {
        return merge_append_values(vcf_val, db_val);
    }

    if let Some(val) = vcf_value {
        if !is_missing_value(val) || mode.carry_over_missing {
            return val.to_string();
        }
    }

    db_val.to_string()
}
