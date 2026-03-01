use anyhow::{Result, anyhow};
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs;

#[derive(Clone, Copy, PartialEq, Eq)]
enum FieldScope {
    Info,
    Format,
    Std,
}

#[derive(Clone, PartialEq, Eq)]
enum IndexSpec {
    All,
    One(usize),
    Range(usize, usize),
    From(usize),
    To(usize),
    List(Vec<usize>),
    Gt,
}

#[derive(Clone)]
struct FieldRef {
    scope: FieldScope,
    name: String,
    sample_sel: Option<IndexSpec>,
    value_sel: Option<IndexSpec>,
}

#[derive(Clone)]
enum Literal {
    Int(i64),
    Float(f64),
    Str(String),
    StrList(HashSet<String>),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BinaryOp {
    Or,
    OrVec,
    And,
    AndVec,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Like,
    NLike,
    Add,
    Sub,
    Mul,
    Div,
    Mod,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum UnaryOp {
    Not,
    Neg,
}

#[derive(Clone, Copy)]
enum FuncName {
    Max,
    Min,
    Avg,
    Sum,
    Median,
    Stdev,
    Abs,
    Count,
    StrLen,
    Phred,
    Binom,
    NPass,
    FPass,
    SMplMax,
    SMplMin,
    SMplAvg,
    SMplSum,
    SMplMedian,
    SMplStdev,
    SMplCount,
}

#[derive(Clone)]
enum Expr {
    Binary {
        op: BinaryOp,
        left: Box<Expr>,
        right: Box<Expr>,
    },
    Unary {
        op: UnaryOp,
        expr: Box<Expr>,
    },
    Literal(Literal),
    Field(FieldRef),
    Func {
        name: FuncName,
        args: Vec<Expr>,
    },
}

#[derive(Clone)]
enum FieldType {
    Integer,
    Float,
    String,
    Flag,
}

#[derive(Clone)]
struct FieldMeta {
    field_type: FieldType,
}

#[derive(Clone)]
struct HeaderMeta {
    info: HashMap<String, FieldMeta>,
    format: HashMap<String, FieldMeta>,
    filters: HashSet<String>,
    samples: Vec<String>,
}

#[derive(Clone)]
enum Value {
    Missing,
    Int(i64),
    Float(f64),
    Str(String),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ValueKind {
    Normal,
    Filter,
    Gt,
    Type,
}

#[derive(Clone)]
struct ValueVec {
    values: Vec<Value>,
    is_str: bool,
    kind: ValueKind,
}

#[derive(Clone)]
struct SampleValues {
    values: Vec<Vec<Value>>,
    mask: Vec<bool>,
    is_str: bool,
    kind: ValueKind,
}

#[derive(Clone)]
enum EvalValue {
    Scalar(ValueVec),
    Samples(SampleValues),
    StrList(HashSet<String>),
}

#[derive(Clone)]
pub struct EvalResult {
    pub pass_site: bool,
    pub pass_samples: Option<Vec<bool>>,
}

pub struct FilterEngine {
    expr: Option<Expr>,
    exclude: bool,
    header: HeaderMeta,
    needed: Option<NeededFields>,
    fast: Option<FastExpr>,
}

impl FilterEngine {
    pub fn new(headers: &[String], expr: Option<&str>, exclude: bool) -> Result<Self> {
        let header = parse_header_meta(headers)?;
        let expr = match expr {
            Some(s) => Some(Parser::new(s, &header).parse_expr()?),
            None => None,
        };
        let needed = expr.as_ref().and_then(collect_needed_fields);
        let fast = expr.as_ref().and_then(|e| build_fast_expr(e, &header));
        Ok(Self {
            expr,
            exclude,
            header,
            needed,
            fast,
        })
    }

    pub fn eval(&self, rec: &crate::vcf::VcfRecord) -> Result<EvalResult> {
        let mut ctx = EvalContext::new(rec, &self.header, self.needed.as_ref());
        let mut result = if let Some(fast) = &self.fast {
            EvalResult {
                pass_site: fast_eval(fast, rec),
                pass_samples: None,
            }
        } else if let Some(expr) = &self.expr {
            eval_bool(expr, &mut ctx)?
        } else {
            EvalResult {
                pass_site: true,
                pass_samples: None,
            }
        };
        if self.exclude {
            result.pass_site = !result.pass_site;
            if let Some(samples) = result.pass_samples.as_mut() {
                for v in samples.iter_mut() {
                    *v = !*v;
                }
            }
        }
        Ok(result)
    }

    pub fn header(&self) -> &HeaderMeta {
        &self.header
    }
}

#[derive(Clone)]
struct NeededFields {
    info_keys: Vec<String>,
    info_index: HashMap<String, usize>,
}

struct EvalContext<'a> {
    rec: &'a crate::vcf::VcfRecord,
    header: &'a HeaderMeta,
    needed: Option<&'a NeededFields>,
    info_cache: Option<HashMap<String, String>>,
    info_values: Option<Vec<Option<String>>>,
    format_cache: Option<FormatCache<'a>>,
}

impl<'a> EvalContext<'a> {
    fn new(
        rec: &'a crate::vcf::VcfRecord,
        header: &'a HeaderMeta,
        needed: Option<&'a NeededFields>,
    ) -> Self {
        Self {
            rec,
            header,
            needed,
            info_cache: None,
            info_values: None,
            format_cache: None,
        }
    }

    fn info_map(&mut self) -> &HashMap<String, String> {
        if self.info_cache.is_none() {
            let map = crate::filter_arch::get_arch().parse_info_simd(&self.rec.info);
            self.info_cache = Some(map);
        }
        self.info_cache.as_ref().unwrap()
    }

    fn info_value(&mut self, key: &str) -> Option<&str> {
        let Some(needed) = self.needed else {
            return self.info_map().get(key).map(|s| s.as_str());
        };
        let Some(&idx) = needed.info_index.get(key) else {
            return None;
        };
        if self.info_values.is_none() {
            let mut values = vec![None; needed.info_keys.len()];
            for item in self.rec.info.split(';') {
                if item.is_empty() {
                    continue;
                }
                if let Some((k, v)) = item.split_once('=') {
                    if let Some(&kidx) = needed.info_index.get(k) {
                        values[kidx] = Some(v.to_string());
                    }
                } else if let Some(&kidx) = needed.info_index.get(item) {
                    values[kidx] = Some(String::new());
                }
            }
            self.info_values = Some(values);
        }
        self.info_values
            .as_ref()
            .and_then(|v| v.get(idx))
            .and_then(|v| v.as_ref().map(|s| s.as_str()))
    }

    fn format_cache(&mut self) -> Option<&FormatCache<'a>> {
        if self.rec.format.is_none() {
            return None;
        }
        if self.format_cache.is_none() {
            let cache = FormatCache::new(self.rec);
            self.format_cache = Some(cache);
        }
        self.format_cache.as_ref()
    }
}

struct FormatCache<'a> {
    keys: Vec<&'a str>,
    key_index: HashMap<&'a str, usize>,
    samples: Vec<Vec<&'a str>>,
}

impl<'a> FormatCache<'a> {
    fn new(rec: &'a crate::vcf::VcfRecord) -> Self {
        let fmt = rec.format.as_deref().unwrap_or(".");
        let keys: Vec<&str> = fmt.split(':').collect();
        let mut key_index = HashMap::new();
        for (i, k) in keys.iter().enumerate() {
            key_index.insert(*k, i);
        }
        let mut samples = Vec::with_capacity(rec.samples.len());
        for s in &rec.samples {
            let vals: Vec<&str> = s.split(':').collect();
            samples.push(vals);
        }
        Self {
            keys,
            key_index,
            samples,
        }
    }

    fn get_value(&self, sample_idx: usize, key: &str) -> Option<&str> {
        let idx = *self.key_index.get(key)?;
        self.samples
            .get(sample_idx)
            .and_then(|v| v.get(idx))
            .copied()
    }

    fn has_key(&self, key: &str) -> bool {
        self.key_index.contains_key(key)
    }
}

fn parse_header_meta(headers: &[String]) -> Result<HeaderMeta> {
    let mut info = HashMap::new();
    let mut format = HashMap::new();
    let mut filters = HashSet::new();
    let mut samples = Vec::new();

    for h in headers {
        if h.starts_with("##INFO=") {
            if let Some((id, ty)) = parse_header_kv(h, "ID", "Type") {
                info.insert(id, FieldMeta { field_type: ty });
            }
        } else if h.starts_with("##FORMAT=") {
            if let Some((id, ty)) = parse_header_kv(h, "ID", "Type") {
                format.insert(id, FieldMeta { field_type: ty });
            }
        } else if h.starts_with("##FILTER=") {
            if let Some(id) = extract_header_id(h, "ID") {
                filters.insert(id);
            }
        } else if h.starts_with("#CHROM") {
            let parts: Vec<&str> = h.split('\t').collect();
            if parts.len() > 9 {
                samples = parts[9..].iter().map(|s| s.to_string()).collect();
            }
        }
    }

    Ok(HeaderMeta {
        info,
        format,
        filters,
        samples,
    })
}

fn parse_header_kv(line: &str, key_id: &str, key_type: &str) -> Option<(String, FieldType)> {
    let id = extract_header_id(line, key_id)?;
    let ty = extract_header_id(line, key_type)?;
    let field_type = match ty.as_str() {
        "Integer" => FieldType::Integer,
        "Float" => FieldType::Float,
        "String" => FieldType::String,
        "Flag" => FieldType::Flag,
        _ => FieldType::String,
    };
    Some((id, field_type))
}

fn extract_header_id(line: &str, key: &str) -> Option<String> {
    let pat = format!("{}=", key);
    let start = line.find(&pat)? + pat.len();
    let rest = &line[start..];
    let end = rest.find(',').or_else(|| rest.find('>'))?;
    Some(rest[..end].to_string())
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    header: &'a HeaderMeta,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str, header: &'a HeaderMeta) -> Self {
        Self {
            chars: input.chars().collect(),
            pos: 0,
            header,
        }
    }

    fn parse_expr(&mut self) -> Result<Expr> {
        let expr = self.parse_or()?;
        self.skip_ws();
        if self.pos < self.chars.len() {
            return Err(anyhow!("Unexpected input at {}", self.pos));
        }
        Ok(expr)
    }

    fn parse_or(&mut self) -> Result<Expr> {
        let mut node = self.parse_and()?;
        loop {
            self.skip_ws();
            if self.consume_str("||") {
                let right = self.parse_and()?;
                node = Expr::Binary {
                    op: BinaryOp::OrVec,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else if self.consume_char('|') {
                let right = self.parse_and()?;
                node = Expr::Binary {
                    op: BinaryOp::Or,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else {
                break;
            }
        }
        Ok(node)
    }

    fn parse_and(&mut self) -> Result<Expr> {
        let mut node = self.parse_cmp()?;
        loop {
            self.skip_ws();
            if self.consume_str("&&") {
                let right = self.parse_cmp()?;
                node = Expr::Binary {
                    op: BinaryOp::AndVec,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else if self.consume_char('&') {
                let right = self.parse_cmp()?;
                node = Expr::Binary {
                    op: BinaryOp::And,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else {
                break;
            }
        }
        Ok(node)
    }

    fn parse_cmp(&mut self) -> Result<Expr> {
        let mut node = self.parse_add()?;
        self.skip_ws();
        if let Some(op) = self.parse_cmp_op() {
            let right = self.parse_add()?;
            node = Expr::Binary {
                op,
                left: Box::new(node),
                right: Box::new(right),
            };
        }
        Ok(node)
    }

    fn parse_add(&mut self) -> Result<Expr> {
        let mut node = self.parse_mul()?;
        loop {
            self.skip_ws();
            if self.consume_char('+') {
                let right = self.parse_mul()?;
                node = Expr::Binary {
                    op: BinaryOp::Add,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else if self.consume_char('-') {
                let right = self.parse_mul()?;
                node = Expr::Binary {
                    op: BinaryOp::Sub,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else {
                break;
            }
        }
        Ok(node)
    }

    fn parse_mul(&mut self) -> Result<Expr> {
        let mut node = self.parse_unary()?;
        loop {
            self.skip_ws();
            if self.consume_char('*') {
                let right = self.parse_unary()?;
                node = Expr::Binary {
                    op: BinaryOp::Mul,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else if self.consume_char('/') {
                let right = self.parse_unary()?;
                node = Expr::Binary {
                    op: BinaryOp::Div,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else if self.consume_char('%') {
                let right = self.parse_unary()?;
                node = Expr::Binary {
                    op: BinaryOp::Mod,
                    left: Box::new(node),
                    right: Box::new(right),
                };
            } else {
                break;
            }
        }
        Ok(node)
    }

    fn parse_unary(&mut self) -> Result<Expr> {
        self.skip_ws();
        if self.consume_char('!') {
            let expr = self.parse_unary()?;
            return Ok(Expr::Unary {
                op: UnaryOp::Not,
                expr: Box::new(expr),
            });
        }
        if self.consume_char('-') {
            let expr = self.parse_unary()?;
            return Ok(Expr::Unary {
                op: UnaryOp::Neg,
                expr: Box::new(expr),
            });
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Expr> {
        self.skip_ws();
        if self.consume_char('(') {
            let expr = self.parse_or()?;
            self.skip_ws();
            if !self.consume_char(')') {
                return Err(anyhow!("Missing )"));
            }
            return Ok(expr);
        }
        if let Some(ch) = self.peek_char() {
            if ch == '"' || ch == '\'' {
                let s = self.parse_string()?;
                return Ok(Expr::Literal(Literal::Str(s)));
            }
            if ch.is_ascii_digit() || (ch == '.' && self.peek_next_is_digit()) {
                let lit = self.parse_number()?;
                return Ok(Expr::Literal(lit));
            }
            if ch == '@' {
                let list = self.parse_file_list()?;
                return Ok(Expr::Literal(Literal::StrList(list)));
            }
            let ident = self.parse_ident()?;
            let ident = self.maybe_qualified_ident(ident)?;
            if self.consume_char('(') {
                let args = self.parse_args()?;
                let name = parse_func_name(&ident)?;
                return Ok(Expr::Func { name, args });
            }
            let field = self.parse_field_from_ident(ident)?;
            return Ok(Expr::Field(field));
        }
        Err(anyhow!("Unexpected end of input"))
    }

    fn parse_args(&mut self) -> Result<Vec<Expr>> {
        let mut args = Vec::new();
        self.skip_ws();
        if self.consume_char(')') {
            return Ok(args);
        }
        loop {
            let expr = self.parse_or()?;
            args.push(expr);
            self.skip_ws();
            if self.consume_char(')') {
                break;
            }
            if !self.consume_char(',') {
                return Err(anyhow!("Expected , or )"));
            }
        }
        Ok(args)
    }

    fn parse_field_from_ident(&mut self, ident: String) -> Result<FieldRef> {
        let scope;
        let mut name = ident;
        let mut sample_sel = None;
        let mut value_sel = None;

        let upper = name.to_ascii_uppercase();
        if upper == "CHROM"
            || upper == "POS"
            || upper == "ID"
            || upper == "REF"
            || upper == "ALT"
            || upper == "QUAL"
            || upper == "FILTER"
            || upper == "TYPE"
            || upper == "N_ALT"
            || upper == "N_SAMPLES"
            || upper == "MAC"
            || upper == "%ILEN"
            || upper == "ILEN"
            || upper == "F_MISSING"
        {
            scope = FieldScope::Std;
            name = upper;
        } else if upper.starts_with("INFO/") {
            scope = FieldScope::Info;
            name = name[5..].to_string();
        } else if upper.starts_with("FORMAT/") {
            scope = FieldScope::Format;
            name = name[7..].to_string();
        } else if upper.starts_with("FMT/") {
            scope = FieldScope::Format;
            name = name[4..].to_string();
        } else if self.header.format.contains_key(&name) && !self.header.info.contains_key(&name) {
            scope = FieldScope::Format;
        } else {
            scope = FieldScope::Info;
        }

        self.skip_ws();
        if self.consume_char('[') {
            let spec = self.read_bracket_content()?;
            let (a, b) = parse_index_pair(&spec)?;
            match scope {
                FieldScope::Format => {
                    if b.is_some() {
                        sample_sel = a;
                        value_sel = b;
                    } else {
                        sample_sel = a;
                    }
                }
                _ => {
                    value_sel = a;
                }
            }
        }

        Ok(FieldRef {
            scope,
            name,
            sample_sel,
            value_sel,
        })
    }

    fn read_bracket_content(&mut self) -> Result<String> {
        let mut out = String::new();
        while let Some(ch) = self.peek_char() {
            if ch == ']' {
                self.pos += 1;
                return Ok(out);
            }
            out.push(ch);
            self.pos += 1;
        }
        Err(anyhow!("Missing ]"))
    }

    fn parse_cmp_op(&mut self) -> Option<BinaryOp> {
        self.skip_ws();
        if self.consume_str("==") {
            return Some(BinaryOp::Eq);
        }
        if self.consume_str("!=") {
            return Some(BinaryOp::Ne);
        }
        if self.consume_str("<=") {
            return Some(BinaryOp::Le);
        }
        if self.consume_str(">=") {
            return Some(BinaryOp::Ge);
        }
        if self.consume_str("!~") {
            return Some(BinaryOp::NLike);
        }
        if self.consume_str("=") {
            return Some(BinaryOp::Eq);
        }
        if self.consume_str("<") {
            return Some(BinaryOp::Lt);
        }
        if self.consume_str(">") {
            return Some(BinaryOp::Gt);
        }
        if self.consume_str("~") {
            return Some(BinaryOp::Like);
        }
        None
    }

    fn parse_string(&mut self) -> Result<String> {
        let quote = self.peek_char().unwrap();
        self.pos += 1;
        let mut out = String::new();
        while let Some(ch) = self.peek_char() {
            self.pos += 1;
            if ch == quote {
                return Ok(out);
            }
            if ch == '\\' {
                if let Some(next) = self.peek_char() {
                    self.pos += 1;
                    out.push(next);
                }
            } else {
                out.push(ch);
            }
        }
        Err(anyhow!("Missing string quote"))
    }

    fn parse_number(&mut self) -> Result<Literal> {
        let start = self.pos;
        let mut has_dot = false;
        let mut has_exp = false;
        while let Some(ch) = self.peek_char() {
            if ch.is_ascii_digit() {
                self.pos += 1;
                continue;
            }
            if ch == '.' && !has_dot && !has_exp {
                has_dot = true;
                self.pos += 1;
                continue;
            }
            if (ch == 'e' || ch == 'E') && !has_exp {
                has_exp = true;
                self.pos += 1;
                if let Some(sign) = self.peek_char() {
                    if sign == '+' || sign == '-' {
                        self.pos += 1;
                    }
                }
                continue;
            }
            break;
        }
        let s: String = self.chars[start..self.pos].iter().collect();
        if !has_dot && !has_exp {
            if let Ok(i) = s.parse::<i64>() {
                return Ok(Literal::Int(i));
            }
        }
        let f = s.parse::<f64>()?;
        Ok(Literal::Float(f))
    }

    fn parse_ident(&mut self) -> Result<String> {
        let start = self.pos;
        while let Some(ch) = self.peek_char() {
            if ch.is_ascii_alphanumeric() || ch == '_' || ch == '.' || ch == '%' {
                self.pos += 1;
                continue;
            }
            break;
        }
        if self.pos == start {
            return Err(anyhow!("Expected identifier"));
        }
        Ok(self.chars[start..self.pos].iter().collect())
    }

    fn maybe_qualified_ident(&mut self, ident: String) -> Result<String> {
        let upper = ident.to_ascii_uppercase();
        if (upper == "INFO" || upper == "FORMAT" || upper == "FMT") && self.consume_char('/') {
            let tail = self.parse_ident()?;
            return Ok(format!("{}/{}", ident, tail));
        }
        Ok(ident)
    }

    fn parse_file_list(&mut self) -> Result<HashSet<String>> {
        self.pos += 1;
        let start = self.pos;
        while let Some(ch) = self.peek_char() {
            if ch.is_whitespace() || ch == ')' || ch == ']' || ch == '=' || ch == '!' {
                break;
            }
            self.pos += 1;
        }
        let path: String = self.chars[start..self.pos].iter().collect();
        let content = fs::read_to_string(path)?;
        let mut set = HashSet::new();
        for line in content.lines() {
            let v = line.trim();
            if !v.is_empty() {
                set.insert(v.to_string());
            }
        }
        Ok(set)
    }

    fn skip_ws(&mut self) {
        while let Some(ch) = self.peek_char() {
            if ch.is_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn peek_char(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn consume_char(&mut self, ch: char) -> bool {
        if self.peek_char() == Some(ch) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn consume_str(&mut self, s: &str) -> bool {
        let len = s.chars().count();
        if self.pos + len > self.chars.len() {
            return false;
        }
        let segment: String = self.chars[self.pos..self.pos + len].iter().collect();
        if segment == s {
            self.pos += len;
            true
        } else {
            false
        }
    }

    fn peek_next_is_digit(&self) -> bool {
        self.chars
            .get(self.pos + 1)
            .map(|c| c.is_ascii_digit())
            .unwrap_or(false)
    }
}

fn parse_func_name(name: &str) -> Result<FuncName> {
    let u = name.to_ascii_uppercase();
    let func = match u.as_str() {
        "MAX" => FuncName::Max,
        "MIN" => FuncName::Min,
        "MEAN" | "AVG" => FuncName::Avg,
        "SUM" => FuncName::Sum,
        "MEDIAN" => FuncName::Median,
        "STDEV" => FuncName::Stdev,
        "ABS" => FuncName::Abs,
        "COUNT" => FuncName::Count,
        "STRLEN" => FuncName::StrLen,
        "PHRED" => FuncName::Phred,
        "BINOM" => FuncName::Binom,
        "N_PASS" => FuncName::NPass,
        "F_PASS" => FuncName::FPass,
        "SMPL_MAX" | "SMPLMAX" | "SMAX" => FuncName::SMplMax,
        "SMPL_MIN" | "SMPLMIN" | "SMIN" => FuncName::SMplMin,
        "SMPL_AVG" | "SMPL_MEAN" | "SMPLAVG" | "SMEAN" | "SAVG" => FuncName::SMplAvg,
        "SMPL_SUM" | "SMPLSUM" | "SSUM" => FuncName::SMplSum,
        "SMPL_MEDIAN" | "SMPLMEDIAN" | "SMEDIAN" => FuncName::SMplMedian,
        "SMPL_STDEV" | "SMPLSTDEV" | "SSTDEV" => FuncName::SMplStdev,
        "SMPL_COUNT" | "SMPLCOUNT" | "SCOUNT" => FuncName::SMplCount,
        _ => return Err(anyhow!("Unknown function {}", name)),
    };
    Ok(func)
}

fn parse_index_pair(spec: &str) -> Result<(Option<IndexSpec>, Option<IndexSpec>)> {
    let mut parts = spec.splitn(2, ':');
    let left = parts.next().unwrap_or("");
    let right = parts.next();
    let a = if left.is_empty() {
        Some(IndexSpec::All)
    } else {
        Some(parse_index_spec(left)?)
    };
    let b = if let Some(r) = right {
        if r.is_empty() {
            Some(IndexSpec::All)
        } else {
            Some(parse_index_spec(r)?)
        }
    } else {
        None
    };
    Ok((a, b))
}

fn parse_index_spec(spec: &str) -> Result<IndexSpec> {
    let s = spec.trim();
    if s == "*" {
        return Ok(IndexSpec::All);
    }
    let upper = s.to_ascii_uppercase();
    if upper == "GT" {
        return Ok(IndexSpec::Gt);
    }
    if let Some((a, b)) = s.split_once("..") {
        if a.is_empty() {
            let end = b.parse::<usize>()?;
            return Ok(IndexSpec::To(end));
        }
        if b.is_empty() {
            let start = a.parse::<usize>()?;
            return Ok(IndexSpec::From(start));
        }
        let start = a.parse::<usize>()?;
        let end = b.parse::<usize>()?;
        return Ok(IndexSpec::Range(start, end));
    }
    if let Some((a, b)) = s.split_once('-') {
        let start = a.parse::<usize>()?;
        let end = b.parse::<usize>()?;
        return Ok(IndexSpec::Range(start, end));
    }
    if s.contains(',') {
        let mut out = Vec::new();
        for part in s.split(',') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            if let Some((a, b)) = part.split_once('-') {
                let start = a.parse::<usize>()?;
                let end = b.parse::<usize>()?;
                for v in start..=end {
                    out.push(v);
                }
            } else {
                out.push(part.parse::<usize>()?);
            }
        }
        return Ok(IndexSpec::List(out));
    }
    let v = s.parse::<usize>()?;
    Ok(IndexSpec::One(v))
}

fn collect_needed_fields(expr: &Expr) -> Option<NeededFields> {
    let mut info = HashSet::new();
    fn walk(expr: &Expr, info: &mut HashSet<String>) {
        match expr {
            Expr::Field(f) => {
                if f.scope == FieldScope::Info {
                    info.insert(f.name.clone());
                }
            }
            Expr::Binary { left, right, .. } => {
                walk(left, info);
                walk(right, info);
            }
            Expr::Unary { expr, .. } => walk(expr, info),
            Expr::Func { args, .. } => {
                for a in args {
                    walk(a, info);
                }
            }
            Expr::Literal(_) => {}
        }
    }
    walk(expr, &mut info);
    if info.is_empty() {
        return None;
    }
    let mut info_keys: Vec<String> = info.into_iter().collect();
    info_keys.sort();
    let mut info_index = HashMap::new();
    for (i, k) in info_keys.iter().enumerate() {
        info_index.insert(k.clone(), i);
    }
    Some(NeededFields {
        info_keys,
        info_index,
    })
}

#[derive(Clone, Copy)]
enum FastField {
    Qual,
    Info(usize),
    InfoFlag(usize),
}

#[derive(Clone, Copy)]
struct FastTerm {
    field: FastField,
    op: BinaryOp,
    value: f64,
}

#[derive(Clone)]
struct FastExpr {
    terms: Vec<FastTerm>,
    info_keys: Vec<String>,
}

fn build_fast_expr(expr: &Expr, header: &HeaderMeta) -> Option<FastExpr> {
    let mut terms = Vec::new();
    let mut info_keys = Vec::new();
    let mut info_index = HashMap::new();
    if !collect_fast_terms(expr, header, &mut terms, &mut info_keys, &mut info_index) {
        return None;
    }
    Some(FastExpr { terms, info_keys })
}

fn collect_fast_terms(
    expr: &Expr,
    header: &HeaderMeta,
    terms: &mut Vec<FastTerm>,
    info_keys: &mut Vec<String>,
    info_index: &mut HashMap<String, usize>,
) -> bool {
    match expr {
        Expr::Binary { op, left, right } => {
            if matches!(op, BinaryOp::And | BinaryOp::AndVec) {
                return collect_fast_terms(left, header, terms, info_keys, info_index)
                    && collect_fast_terms(right, header, terms, info_keys, info_index);
            }
            if matches!(
                op,
                BinaryOp::Eq
                    | BinaryOp::Ne
                    | BinaryOp::Lt
                    | BinaryOp::Le
                    | BinaryOp::Gt
                    | BinaryOp::Ge
            ) {
                return collect_fast_term(op, left, right, header, terms, info_keys, info_index)
                    || collect_fast_term(op, right, left, header, terms, info_keys, info_index);
            }
            false
        }
        _ => false,
    }
}

fn collect_fast_term(
    op: &BinaryOp,
    field_expr: &Expr,
    lit_expr: &Expr,
    header: &HeaderMeta,
    terms: &mut Vec<FastTerm>,
    info_keys: &mut Vec<String>,
    info_index: &mut HashMap<String, usize>,
) -> bool {
    let lit = match lit_expr {
        Expr::Literal(Literal::Int(i)) => *i as f64,
        Expr::Literal(Literal::Float(f)) => *f,
        _ => return false,
    };
    let field = match field_expr {
        Expr::Field(f) => f,
        _ => return false,
    };
    if field.sample_sel.is_some() || field.value_sel.is_some() {
        return false;
    }
    let term_field = match field.scope {
        FieldScope::Std => {
            if field.name == "QUAL" {
                FastField::Qual
            } else {
                return false;
            }
        }
        FieldScope::Info => {
            let meta = header.info.get(&field.name);
            let field_type = meta
                .map(|m| m.field_type.clone())
                .unwrap_or(FieldType::String);
            match field_type {
                FieldType::Integer | FieldType::Float => {
                    let idx = *info_index.entry(field.name.clone()).or_insert_with(|| {
                        let i = info_keys.len();
                        info_keys.push(field.name.clone());
                        i
                    });
                    FastField::Info(idx)
                }
                FieldType::Flag => {
                    let idx = *info_index.entry(field.name.clone()).or_insert_with(|| {
                        let i = info_keys.len();
                        info_keys.push(field.name.clone());
                        i
                    });
                    FastField::InfoFlag(idx)
                }
                FieldType::String => return false,
            }
        }
        FieldScope::Format => return false,
    };
    terms.push(FastTerm {
        field: term_field,
        op: *op,
        value: lit,
    });
    true
}

fn fast_eval(fast: &FastExpr, rec: &crate::vcf::VcfRecord) -> bool {
    for term in &fast.terms {
        let pass = match term.field {
            FastField::Qual => fast_cmp_values(term.op, parse_qual(&rec.qual), term.value),
            FastField::Info(idx) => {
                let key = &fast.info_keys[idx];
                fast_cmp_values(term.op, info_numbers(rec.info.as_str(), key), term.value)
            }
            FastField::InfoFlag(idx) => {
                let key = &fast.info_keys[idx];
                let present = info_has_flag(&rec.info, key);
                let values = if present { vec![1.0] } else { Vec::new() };
                fast_cmp_values(term.op, values, term.value)
            }
        };
        if !pass {
            return false;
        }
    }
    true
}

fn parse_qual(qual: &str) -> Vec<f64> {
    if qual.is_empty() || qual == "." {
        return Vec::new();
    }
    if let Ok(v) = qual.parse::<f64>() {
        vec![v]
    } else {
        Vec::new()
    }
}

fn info_numbers(info: &str, key: &str) -> Vec<f64> {
    if info.is_empty() || info == "." {
        return Vec::new();
    }
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        if let Some((k, v)) = item.split_once('=') {
            if k == key {
                return parse_number_list(v);
            }
        }
    }
    Vec::new()
}

fn parse_number_list(s: &str) -> Vec<f64> {
    let mut out = Vec::new();
    for part in s.split(',') {
        if part == "." || part.is_empty() {
            out.push(f64::NAN);
            continue;
        }
        if let Ok(v) = part.parse::<f64>() {
            out.push(v);
        }
    }
    out
}

fn fast_cmp_values(op: BinaryOp, values: Vec<f64>, rhs: f64) -> bool {
    let miss_one = matches!(op, BinaryOp::Ne);
    if values.is_empty() {
        return miss_one;
    }
    let mut any_value = false;
    for v in values {
        if v.is_nan() {
            if miss_one {
                return true;
            }
            continue;
        }
        any_value = true;
        let ok = match op {
            BinaryOp::Eq => v == rhs,
            BinaryOp::Ne => v != rhs,
            BinaryOp::Lt => v < rhs,
            BinaryOp::Le => v <= rhs,
            BinaryOp::Gt => v > rhs,
            BinaryOp::Ge => v >= rhs,
            _ => false,
        };
        if ok {
            return true;
        }
    }
    if any_value { false } else { miss_one }
}
fn eval_bool(expr: &Expr, ctx: &mut EvalContext) -> Result<EvalResult> {
    match expr {
        Expr::Binary { op, left, right } => match op {
            BinaryOp::Or | BinaryOp::OrVec | BinaryOp::And | BinaryOp::AndVec => {
                let a = eval_bool(left, ctx)?;
                let b = eval_bool(right, ctx)?;
                Ok(combine_logic(*op, a, b))
            }
            BinaryOp::Eq
            | BinaryOp::Ne
            | BinaryOp::Lt
            | BinaryOp::Le
            | BinaryOp::Gt
            | BinaryOp::Ge
            | BinaryOp::Like
            | BinaryOp::NLike => {
                let a = eval_value(left, ctx)?;
                let b = eval_value(right, ctx)?;
                Ok(compare_values(*op, a, b, ctx)?)
            }
            BinaryOp::Add | BinaryOp::Sub | BinaryOp::Mul | BinaryOp::Div | BinaryOp::Mod => {
                let v = eval_value(expr, ctx)?;
                Ok(value_to_bool(v))
            }
        },
        Expr::Unary { op, expr } => {
            let mut v = eval_bool(expr, ctx)?;
            match op {
                UnaryOp::Not => {
                    v.pass_site = !v.pass_site;
                    if let Some(samples) = v.pass_samples.as_mut() {
                        for s in samples.iter_mut() {
                            *s = !*s;
                        }
                    }
                }
                UnaryOp::Neg => {}
            }
            Ok(v)
        }
        Expr::Literal(_) | Expr::Field(_) | Expr::Func { .. } => {
            let v = eval_value(expr, ctx)?;
            Ok(value_to_bool(v))
        }
    }
}

fn eval_value(expr: &Expr, ctx: &mut EvalContext) -> Result<EvalValue> {
    match expr {
        Expr::Literal(lit) => match lit {
            Literal::Int(i) => Ok(EvalValue::Scalar(ValueVec {
                values: vec![Value::Int(*i)],
                is_str: false,
                kind: ValueKind::Normal,
            })),
            Literal::Float(f) => Ok(EvalValue::Scalar(ValueVec {
                values: vec![Value::Float(*f)],
                is_str: false,
                kind: ValueKind::Normal,
            })),
            Literal::Str(s) => {
                if s == "." {
                    Ok(EvalValue::Scalar(ValueVec {
                        values: vec![Value::Missing],
                        is_str: true,
                        kind: ValueKind::Normal,
                    }))
                } else {
                    Ok(EvalValue::Scalar(ValueVec {
                        values: vec![Value::Str(s.clone())],
                        is_str: true,
                        kind: ValueKind::Normal,
                    }))
                }
            }
            Literal::StrList(s) => Ok(EvalValue::StrList(s.clone())),
        },
        Expr::Field(f) => eval_field(f, ctx),
        Expr::Unary { op, expr } => {
            let v = eval_value(expr, ctx)?;
            match op {
                UnaryOp::Neg => Ok(apply_unary_neg(v)),
                UnaryOp::Not => Ok(v),
            }
        }
        Expr::Binary { op, left, right } => {
            let a = eval_value(left, ctx)?;
            let b = eval_value(right, ctx)?;
            match op {
                BinaryOp::Add | BinaryOp::Sub | BinaryOp::Mul | BinaryOp::Div | BinaryOp::Mod => {
                    Ok(apply_arith(*op, a, b))
                }
                _ => Ok(value_from_bool(eval_bool(expr, ctx)?)),
            }
        }
        Expr::Func { name, args } => eval_func(*name, args, ctx),
    }
}

fn eval_field(field: &FieldRef, ctx: &mut EvalContext) -> Result<EvalValue> {
    match field.scope {
        FieldScope::Std => eval_std_field(field, ctx),
        FieldScope::Info => eval_info_field(field, ctx),
        FieldScope::Format => eval_format_field(field, ctx),
    }
}

fn eval_std_field(field: &FieldRef, ctx: &mut EvalContext) -> Result<EvalValue> {
    let name = field.name.as_str();
    let rec = ctx.rec;
    match name {
        "CHROM" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Str(rec.chrom.clone())],
            is_str: true,
            kind: ValueKind::Normal,
        })),
        "POS" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Int(rec.pos as i64)],
            is_str: false,
            kind: ValueKind::Normal,
        })),
        "ID" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Str(rec.id.clone())],
            is_str: true,
            kind: ValueKind::Normal,
        })),
        "REF" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Str(rec.ref_allele.clone())],
            is_str: true,
            kind: ValueKind::Normal,
        })),
        "ALT" => Ok(EvalValue::Scalar(ValueVec {
            values: if rec.alt == "." {
                vec![Value::Str(".".to_string())]
            } else {
                split_to_values(&rec.alt, FieldType::String)
            },
            is_str: true,
            kind: ValueKind::Normal,
        })),
        "QUAL" => Ok(EvalValue::Scalar(ValueVec {
            values: parse_numeric_value(&rec.qual),
            is_str: false,
            kind: ValueKind::Normal,
        })),
        "FILTER" => Ok(EvalValue::Scalar(ValueVec {
            values: split_filter_values(&rec.filter),
            is_str: true,
            kind: ValueKind::Filter,
        })),
        "TYPE" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Int(variant_type_mask(rec) as i64)],
            is_str: false,
            kind: ValueKind::Type,
        })),
        "N_ALT" => {
            let n = if rec.alt == "." || rec.alt.is_empty() {
                0
            } else {
                rec.alt
                    .split(',')
                    .filter(|a| !a.is_empty() && *a != ".")
                    .count() as i64
            };
            Ok(EvalValue::Scalar(ValueVec {
                values: vec![Value::Int(n)],
                is_str: false,
                kind: ValueKind::Normal,
            }))
        }
        "N_SAMPLES" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Int(ctx.header.samples.len() as i64)],
            is_str: false,
            kind: ValueKind::Normal,
        })),
        "MAC" => {
            let mut values = Vec::new();
            let ac_s = ctx.info_value("AC").map(|s| s.to_string());
            let an_s = ctx.info_value("AN").map(|s| s.to_string());
            if let Some(ac) = ac_s {
                let an = ctx
                    .info_value("AN")
                    .and_then(|v| v.split(',').next())
                    .and_then(|v| v.trim().parse::<i64>().ok());
                let an = an.or_else(|| {
                    an_s.as_deref()
                        .and_then(|v| v.split(',').next())
                        .and_then(|v| v.trim().parse::<i64>().ok())
                });
                for v in ac.split(',') {
                    if let Ok(ai) = v.trim().parse::<i64>() {
                        let mac = an.map(|x| ai.min(x - ai)).unwrap_or(ai);
                        values.push(Value::Int(mac));
                    }
                }
            }
            if let Some(sel) = &field.value_sel {
                values = select_values(&values, sel);
            }
            Ok(EvalValue::Scalar(ValueVec {
                values,
                is_str: false,
                kind: ValueKind::Normal,
            }))
        }
        "%ILEN" | "ILEN" => Ok(EvalValue::Scalar(ValueVec {
            values: if let Some(v) = ctx
                .info_value("ILEN")
                .and_then(|x| x.split(',').next())
                .and_then(|x| x.trim().parse::<i64>().ok())
            {
                vec![Value::Int(v)]
            } else {
                vec![Value::Int(variant_ilen(rec) as i64)]
            },
            is_str: false,
            kind: ValueKind::Normal,
        })),
        "F_MISSING" => Ok(EvalValue::Scalar(ValueVec {
            values: vec![Value::Float(f_missing(ctx))],
            is_str: false,
            kind: ValueKind::Normal,
        })),
        _ => Ok(EvalValue::Scalar(ValueVec {
            values: Vec::new(),
            is_str: false,
            kind: ValueKind::Normal,
        })),
    }
}

fn eval_info_field(field: &FieldRef, ctx: &mut EvalContext) -> Result<EvalValue> {
    let field_type = ctx
        .header
        .info
        .get(&field.name)
        .map(|m| m.field_type.clone())
        .unwrap_or(FieldType::String);
    let val = ctx.info_value(&field.name);
    let mut values = match val {
        Some(v) => {
            if matches!(field_type, FieldType::Flag) && v.is_empty() {
                vec![Value::Int(1)]
            } else {
                split_to_values(v, field_type.clone())
            }
        }
        None => {
            if matches!(field_type, FieldType::Flag) && info_has_flag(&ctx.rec.info, &field.name) {
                vec![Value::Int(1)]
            } else {
                Vec::new()
            }
        }
    };
    if let Some(sel) = &field.value_sel {
        values = select_values(&values, sel);
    }
    Ok(EvalValue::Scalar(ValueVec {
        values,
        is_str: matches!(field_type, FieldType::String),
        kind: ValueKind::Normal,
    }))
}

fn eval_format_field(field: &FieldRef, ctx: &mut EvalContext) -> Result<EvalValue> {
    let field_type = ctx
        .header
        .format
        .get(&field.name)
        .map(|m| m.field_type.clone())
        .unwrap_or(FieldType::String);
    let nsamples = ctx.header.samples.len();
    let is_gt = field.name == "GT";
    let cache = match ctx.format_cache() {
        Some(c) => c,
        None => {
            return Ok(EvalValue::Samples(SampleValues {
                values: Vec::new(),
                mask: Vec::new(),
                is_str: false,
                kind: ValueKind::Normal,
            }));
        }
    };
    if !cache.has_key(&field.name) {
        return Ok(EvalValue::Samples(SampleValues {
            values: vec![Vec::new(); nsamples],
            mask: vec![true; nsamples],
            is_str: false,
            kind: if is_gt {
                ValueKind::Gt
            } else {
                ValueKind::Normal
            },
        }));
    }
    let mut values = vec![Vec::new(); nsamples];
    let mut mask = vec![true; nsamples];

    let sample_sel = field.sample_sel.clone().unwrap_or(IndexSpec::All);
    let sample_idx = select_sample_indices(nsamples, &sample_sel);
    for i in 0..nsamples {
        if !sample_idx.contains(&i) {
            mask[i] = false;
        }
    }

    for &i in &sample_idx {
        if let Some(raw) = cache.get_value(i, &field.name) {
            let mut vals = if field.name == "GT" {
                parse_gt_value(raw)
            } else {
                split_to_values(raw, field_type.clone())
            };
            if let Some(sel) = &field.value_sel {
                if matches!(sel, IndexSpec::Gt) {
                    vals = select_by_gt(cache, i, &vals);
                } else {
                    vals = select_values(&vals, sel);
                }
            }
            values[i] = vals;
        } else {
            values[i] = vec![Value::Missing];
        }
    }

    Ok(EvalValue::Samples(SampleValues {
        values,
        mask,
        is_str: matches!(field_type, FieldType::String),
        kind: if is_gt {
            ValueKind::Gt
        } else {
            ValueKind::Normal
        },
    }))
}

fn split_filter_values(input: &str) -> Vec<Value> {
    if input.is_empty() || input == "." || input == "PASS" {
        return Vec::new();
    }
    input
        .split(';')
        .filter(|s| !s.is_empty())
        .map(|s| Value::Str(s.to_string()))
        .collect()
}

fn info_has_flag(info: &str, key: &str) -> bool {
    if info.is_empty() || info == "." {
        return false;
    }
    for part in info.split(';') {
        if part == key {
            return true;
        }
    }
    false
}

fn split_to_values(input: &str, field_type: FieldType) -> Vec<Value> {
    if input.is_empty() {
        return Vec::new();
    }
    if input == "." {
        return vec![Value::Missing];
    }
    let parts: Vec<&str> = input.split(',').collect();
    let mut out = Vec::with_capacity(parts.len());
    for p in parts {
        if p == "." {
            out.push(Value::Missing);
            continue;
        }
        match field_type {
            FieldType::Integer => {
                if let Ok(i) = p.parse::<i64>() {
                    out.push(Value::Int(i));
                } else {
                    out.push(Value::Missing);
                }
            }
            FieldType::Float => {
                if let Ok(f) = p.parse::<f64>() {
                    out.push(Value::Float(f));
                } else {
                    out.push(Value::Missing);
                }
            }
            FieldType::String | FieldType::Flag => out.push(Value::Str(p.to_string())),
        }
    }
    out
}

fn parse_numeric_value(input: &str) -> Vec<Value> {
    if input.is_empty() || input == "." {
        return Vec::new();
    }
    if let Ok(v) = input.parse::<f64>() {
        return vec![Value::Float(v)];
    }
    Vec::new()
}

fn parse_gt_value(raw: &str) -> Vec<Value> {
    if raw.is_empty() || raw == "." || raw == "./." || raw == ".|." {
        return vec![Value::Missing];
    }
    vec![Value::Str(raw.to_string())]
}

fn select_values(values: &[Value], sel: &IndexSpec) -> Vec<Value> {
    if values.is_empty() {
        return Vec::new();
    }
    match sel {
        IndexSpec::All => values.to_vec(),
        IndexSpec::One(i) => values.get(*i).map(|v| vec![v.clone()]).unwrap_or_default(),
        IndexSpec::Range(a, b) => values
            .iter()
            .enumerate()
            .filter_map(|(i, v)| {
                if i >= *a && i <= *b {
                    Some(v.clone())
                } else {
                    None
                }
            })
            .collect(),
        IndexSpec::From(a) => values
            .iter()
            .enumerate()
            .filter_map(|(i, v)| if i >= *a { Some(v.clone()) } else { None })
            .collect(),
        IndexSpec::To(b) => values
            .iter()
            .enumerate()
            .filter_map(|(i, v)| if i <= *b { Some(v.clone()) } else { None })
            .collect(),
        IndexSpec::List(list) => list
            .iter()
            .filter_map(|i| values.get(*i).cloned())
            .collect(),
        IndexSpec::Gt => values.to_vec(),
    }
}

fn select_sample_indices(nsamples: usize, sel: &IndexSpec) -> Vec<usize> {
    match sel {
        IndexSpec::All => (0..nsamples).collect(),
        IndexSpec::One(i) => {
            if *i < nsamples {
                vec![*i]
            } else {
                Vec::new()
            }
        }
        IndexSpec::Range(a, b) => (*a..=(*b)).filter(|i| *i < nsamples).collect(),
        IndexSpec::From(a) => (*a..nsamples).collect(),
        IndexSpec::To(b) => (0..=(*b)).filter(|i| *i < nsamples).collect(),
        IndexSpec::List(list) => list.iter().copied().filter(|i| *i < nsamples).collect(),
        IndexSpec::Gt => (0..nsamples).collect(),
    }
}

fn select_by_gt(cache: &FormatCache, sample_idx: usize, values: &[Value]) -> Vec<Value> {
    let gt = parse_gt_from_cache(cache, sample_idx);
    if gt.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    for a in gt {
        if a >= 0 {
            let idx = a as usize;
            if let Some(v) = values.get(idx) {
                out.push(v.clone());
            }
        }
    }
    out
}

fn parse_gt_from_cache(cache: &FormatCache, sample_idx: usize) -> Vec<i32> {
    let raw = match cache.get_value(sample_idx, "GT") {
        Some(v) => v,
        None => return Vec::new(),
    };
    let mut out = Vec::new();
    for part in raw.split(|c| c == '/' || c == '|') {
        if part == "." || part.is_empty() {
            out.push(-1);
        } else if let Ok(v) = part.parse::<i32>() {
            out.push(v);
        }
    }
    out
}

fn parse_gt(ctx: &mut EvalContext, sample_idx: usize) -> Vec<i32> {
    let cache = match ctx.format_cache() {
        Some(c) => c,
        None => return Vec::new(),
    };
    parse_gt_from_cache(cache, sample_idx)
}

const VCF_SNP: u32 = 1;
const VCF_MNP: u32 = 2;
const VCF_INDEL: u32 = 4;
const VCF_OTHER: u32 = 8;
const VCF_BND: u32 = 16;
const VCF_OVERLAP: u32 = 32;

fn variant_type_mask(rec: &crate::vcf::VcfRecord) -> u32 {
    let ref_bytes = rec.ref_allele.as_bytes();
    let mut mask = 0u32;
    for alt in rec.alt.split(',') {
        let alt = alt.trim();
        if alt.is_empty() {
            continue;
        }
        let t = variant_type_for_alt(ref_bytes, alt.as_bytes());
        mask |= t;
    }
    if mask == 0 { 1 } else { mask << 1 }
}

fn variant_type_for_alt(ref_bytes: &[u8], alt_bytes: &[u8]) -> u32 {
    if alt_bytes == b"*" {
        return VCF_OVERLAP;
    }
    if ref_bytes.len() == 1 && alt_bytes.len() == 1 {
        let a = alt_bytes[0];
        let r = ref_bytes[0];
        if a == b'.' || eq_icase(a, r) {
            return 0;
        }
        if a == b'X' || a == b'x' {
            return 0;
        }
        return VCF_SNP;
    }
    if alt_bytes.first() == Some(&b'<') {
        if alt_bytes.len() >= 3 {
            if (alt_bytes[1] == b'X' || alt_bytes[1] == b'x') && alt_bytes[2] == b'>' {
                return 0;
            }
            if alt_bytes[1] == b'*' && alt_bytes[2] == b'>' {
                return 0;
            }
            if alt_bytes.len() >= 9 && alt_bytes[1..].starts_with(b"NON_REF>") {
                return 0;
            }
        }
        return VCF_OTHER;
    }
    if alt_bytes[0] == b']' || alt_bytes[0] == b'[' {
        return VCF_BND;
    }

    let mut r_i = 0usize;
    let mut a_i = 0usize;
    while r_i < ref_bytes.len() && a_i < alt_bytes.len() && eq_icase(ref_bytes[r_i], alt_bytes[a_i])
    {
        r_i += 1;
        a_i += 1;
    }
    if a_i < alt_bytes.len() && r_i == ref_bytes.len() {
        if alt_bytes[alt_bytes.len() - 1] == b']' || alt_bytes[alt_bytes.len() - 1] == b'[' {
            return VCF_BND;
        }
        return VCF_INDEL;
    } else if r_i < ref_bytes.len() && a_i == alt_bytes.len() {
        return VCF_INDEL;
    } else if r_i == ref_bytes.len() && a_i == alt_bytes.len() {
        return 0;
    }

    let mut re = ref_bytes.len() - 1;
    let mut ae = alt_bytes.len() - 1;
    if alt_bytes[ae] == b']' || alt_bytes[ae] == b'[' {
        return VCF_BND;
    }
    while re > r_i && ae > a_i && eq_icase(ref_bytes[re], alt_bytes[ae]) {
        re -= 1;
        ae -= 1;
    }
    if ae == a_i {
        if re == r_i {
            return VCF_SNP;
        }
        if eq_icase(ref_bytes[re], alt_bytes[ae]) {
            return VCF_INDEL;
        }
        return VCF_OTHER;
    }
    if re == r_i {
        if eq_icase(ref_bytes[re], alt_bytes[ae]) {
            return VCF_INDEL;
        }
        return VCF_OTHER;
    }
    if (re as isize - r_i as isize) == (ae as isize - a_i as isize) {
        VCF_MNP
    } else {
        VCF_OTHER
    }
}

fn eq_icase(a: u8, b: u8) -> bool {
    a.to_ascii_uppercase() == b.to_ascii_uppercase()
}

fn variant_ilen(rec: &crate::vcf::VcfRecord) -> i32 {
    let alts: Vec<&str> = rec.alt.split(',').collect();
    if alts.is_empty() {
        return 0;
    }
    let alt = alts[0];
    (alt.len() as i32) - (rec.ref_allele.len() as i32)
}

fn f_missing(ctx: &mut EvalContext) -> f64 {
    let cache = match ctx.format_cache() {
        Some(c) => c,
        None => return 0.0,
    };
    if !cache.has_key("GT") {
        return 0.0;
    }
    let mut missing = 0;
    for i in 0..ctx.header.samples.len() {
        let gt = parse_gt(ctx, i);
        if gt.iter().all(|v| *v < 0) {
            missing += 1;
        }
    }
    if ctx.header.samples.is_empty() {
        0.0
    } else {
        missing as f64 / ctx.header.samples.len() as f64
    }
}
fn combine_logic(op: BinaryOp, a: EvalResult, b: EvalResult) -> EvalResult {
    match op {
        BinaryOp::Or | BinaryOp::OrVec => logic_or(op, a, b),
        BinaryOp::And | BinaryOp::AndVec => logic_and(op, a, b),
        _ => EvalResult {
            pass_site: false,
            pass_samples: None,
        },
    }
}

fn logic_or(op: BinaryOp, a: EvalResult, b: EvalResult) -> EvalResult {
    if !a.pass_site && !b.pass_site {
        return EvalResult {
            pass_site: false,
            pass_samples: None,
        };
    }
    let mut result = EvalResult {
        pass_site: true,
        pass_samples: None,
    };
    match (a.pass_samples, b.pass_samples) {
        (None, None) => result,
        (Some(sa), None) => {
            if op == BinaryOp::OrVec && !b.pass_site {
                result.pass_samples = Some(sa);
            } else if op == BinaryOp::OrVec {
                result.pass_samples = Some(vec![true; sa.len()]);
            } else {
                result.pass_samples = Some(sa);
            }
            result
        }
        (None, Some(sb)) => {
            if op == BinaryOp::OrVec && !a.pass_site {
                result.pass_samples = Some(sb);
            } else if op == BinaryOp::OrVec {
                result.pass_samples = Some(vec![true; sb.len()]);
            } else {
                result.pass_samples = Some(sb);
            }
            result
        }
        (Some(sa), Some(sb)) => {
            if op == BinaryOp::OrVec {
                let len = sa.len().max(sb.len());
                result.pass_samples = Some(vec![true; len]);
            } else {
                let len = sa.len().min(sb.len());
                let mut out = vec![false; len];
                for i in 0..len {
                    out[i] = sa[i] || sb[i];
                }
                result.pass_samples = Some(out);
            }
            result
        }
    }
}

fn logic_and(op: BinaryOp, a: EvalResult, b: EvalResult) -> EvalResult {
    if !a.pass_site || !b.pass_site {
        return EvalResult {
            pass_site: false,
            pass_samples: None,
        };
    }
    let mut result = EvalResult {
        pass_site: true,
        pass_samples: None,
    };
    match (a.pass_samples, b.pass_samples) {
        (None, None) => result,
        (Some(sa), None) => {
            result.pass_samples = Some(sa);
            result
        }
        (None, Some(sb)) => {
            result.pass_samples = Some(sb);
            result
        }
        (Some(sa), Some(sb)) => {
            let len = sa.len().min(sb.len());
            let mut out = vec![false; len];
            if op == BinaryOp::AndVec {
                for i in 0..len {
                    out[i] = sa[i] || sb[i];
                }
                result.pass_site = true;
            } else {
                let mut any = false;
                for i in 0..len {
                    out[i] = sa[i] && sb[i];
                    if out[i] {
                        any = true;
                    }
                }
                result.pass_site = any;
            }
            result.pass_samples = Some(out);
            result
        }
    }
}

fn compare_values(
    op: BinaryOp,
    a: EvalValue,
    b: EvalValue,
    _ctx: &mut EvalContext,
) -> Result<EvalResult> {
    if let EvalValue::StrList(list) = &a {
        return compare_with_list(op, list, b);
    }
    if let EvalValue::StrList(list) = &b {
        return compare_with_list(op, list, a);
    }
    if let EvalValue::Samples(av) = &a {
        if av.kind == ValueKind::Gt {
            return Ok(compare_gt(op, av, b));
        }
    }
    if let EvalValue::Samples(bv) = &b {
        if bv.kind == ValueKind::Gt {
            return Ok(compare_gt(op, bv, a));
        }
    }
    if let EvalValue::Scalar(av) = &a {
        if av.kind == ValueKind::Type {
            return Ok(compare_type(op, av, b));
        }
    }
    if let EvalValue::Scalar(bv) = &b {
        if bv.kind == ValueKind::Type {
            return Ok(compare_type(op, bv, a));
        }
    }
    if let EvalValue::Scalar(av) = &a {
        if av.kind == ValueKind::Filter {
            return Ok(compare_filter(op, av, b));
        }
    }
    if let EvalValue::Scalar(bv) = &b {
        if bv.kind == ValueKind::Filter {
            return Ok(compare_filter(op, bv, a));
        }
    }
    match (a, b) {
        (EvalValue::Scalar(av), EvalValue::Scalar(bv)) => {
            Ok(compare_scalar(op, &av.values, &bv.values))
        }
        (EvalValue::Samples(av), EvalValue::Scalar(bv)) => {
            Ok(compare_samples_scalar(op, &av, &bv.values))
        }
        (EvalValue::Scalar(av), EvalValue::Samples(bv)) => {
            Ok(compare_samples_scalar(op, &bv, &av.values))
        }
        (EvalValue::Samples(av), EvalValue::Samples(bv)) => Ok(compare_samples(op, &av, &bv)),
        _ => Ok(EvalResult {
            pass_site: false,
            pass_samples: None,
        }),
    }
}

fn compare_with_list(op: BinaryOp, list: &HashSet<String>, other: EvalValue) -> Result<EvalResult> {
    let check = |s: &str| match op {
        BinaryOp::Eq => list.contains(s),
        BinaryOp::Ne => !list.contains(s),
        _ => false,
    };
    match other {
        EvalValue::Scalar(v) => {
            let mut pass = false;
            for val in v.values {
                if let Value::Str(s) = val {
                    if check(&s) {
                        pass = true;
                        break;
                    }
                }
            }
            Ok(EvalResult {
                pass_site: pass,
                pass_samples: None,
            })
        }
        EvalValue::Samples(v) => {
            let mut pass_samples = vec![false; v.values.len()];
            let mut pass_site = false;
            for (i, vals) in v.values.iter().enumerate() {
                if !v.mask[i] {
                    continue;
                }
                for val in vals {
                    if let Value::Str(s) = val {
                        if check(s) {
                            pass_samples[i] = true;
                            pass_site = true;
                            break;
                        }
                    }
                }
            }
            Ok(EvalResult {
                pass_site,
                pass_samples: Some(pass_samples),
            })
        }
        EvalValue::StrList(_) => Ok(EvalResult {
            pass_site: false,
            pass_samples: None,
        }),
    }
}

fn compare_filter(op: BinaryOp, filter: &ValueVec, other: EvalValue) -> EvalResult {
    let mut query = Vec::new();
    match other {
        EvalValue::Scalar(v) => {
            for val in v.values {
                if let Value::Str(s) = val {
                    for part in s.split(';') {
                        if part == "." {
                            continue;
                        }
                        query.push(part.to_string());
                    }
                }
            }
        }
        EvalValue::StrList(list) => {
            for v in list {
                query.push(v);
            }
        }
        _ => {}
    }
    let filter_vals: Vec<String> = filter
        .values
        .iter()
        .filter_map(|v| match v {
            Value::Str(s) => Some(s.clone()),
            _ => None,
        })
        .collect();

    let pass = match op {
        BinaryOp::Like => filter_in(query, filter_vals, true),
        BinaryOp::NLike => filter_in(query, filter_vals, false),
        BinaryOp::Eq => filter_eq(query, filter_vals),
        BinaryOp::Ne => !filter_eq(query, filter_vals),
        _ => false,
    };
    EvalResult {
        pass_site: pass,
        pass_samples: None,
    }
}

fn filter_eq(query: Vec<String>, filter: Vec<String>) -> bool {
    if query.is_empty() && filter.is_empty() {
        return true;
    }
    if query.len() != filter.len() {
        return false;
    }
    for q in &query {
        if !filter.contains(q) {
            return false;
        }
    }
    true
}

fn filter_in(query: Vec<String>, filter: Vec<String>, want_in: bool) -> bool {
    if query.is_empty() {
        return if want_in {
            filter.is_empty()
        } else {
            !filter.is_empty()
        };
    }
    if filter.is_empty() {
        return false;
    }
    if want_in {
        for q in &query {
            if !filter.contains(q) {
                return false;
            }
        }
        true
    } else {
        for q in &query {
            if !filter.contains(q) {
                return true;
            }
        }
        false
    }
}

fn compare_gt(op: BinaryOp, gt: &SampleValues, other: EvalValue) -> EvalResult {
    let query = match other {
        EvalValue::Scalar(v) => extract_gt_query(&v.values),
        EvalValue::StrList(list) => list.into_iter().collect::<Vec<_>>(),
        _ => Vec::new(),
    };
    let mut pass_samples = vec![false; gt.values.len()];
    let mut pass_site = false;
    for i in 0..gt.values.len() {
        if !gt.mask[i] {
            continue;
        }
        let vals = &gt.values[i];
        if vals.is_empty() {
            continue;
        }
        let gt_val = match &vals[0] {
            Value::Str(s) => s.as_str(),
            Value::Missing => ".",
            _ => "",
        };
        let hit = match op {
            BinaryOp::Eq => gt_match_eq(gt_val, &query),
            BinaryOp::Ne => !gt_match_eq(gt_val, &query),
            BinaryOp::Like => gt_match_regex(gt_val, &query, false),
            BinaryOp::NLike => gt_match_regex(gt_val, &query, true),
            _ => false,
        };
        if hit {
            pass_samples[i] = true;
            pass_site = true;
        }
    }
    EvalResult {
        pass_site,
        pass_samples: Some(pass_samples),
    }
}

fn compare_type(op: BinaryOp, type_val: &ValueVec, other: EvalValue) -> EvalResult {
    let tmask = match type_val.values.get(0) {
        Some(Value::Int(i)) => *i as u32,
        Some(Value::Float(f)) => *f as u32,
        _ => 0,
    };
    let mut queries = Vec::new();
    match other {
        EvalValue::Scalar(v) => {
            for val in v.values {
                if let Some(q) = type_mask_from_value(val) {
                    queries.push(q);
                }
            }
        }
        _ => {}
    }
    if queries.is_empty() {
        return EvalResult {
            pass_site: false,
            pass_samples: None,
        };
    }
    let mut pass = false;
    for q in queries {
        let ok = match op {
            BinaryOp::Eq => tmask == q,
            BinaryOp::Ne => tmask != q,
            BinaryOp::Like => (tmask & q) != 0,
            BinaryOp::NLike => (tmask & q) == 0,
            _ => false,
        };
        if ok {
            pass = true;
            break;
        }
    }
    EvalResult {
        pass_site: pass,
        pass_samples: None,
    }
}

fn type_mask_from_value(v: Value) -> Option<u32> {
    match v {
        Value::Int(i) => Some(i as u32),
        Value::Float(f) => Some(f as u32),
        Value::Str(s) => type_mask_from_query(&s),
        Value::Missing => None,
    }
}

fn type_mask_from_query(s: &str) -> Option<u32> {
    let u = s.to_ascii_lowercase();
    match u.as_str() {
        "snp" | "snps" => Some(VCF_SNP << 1),
        "indel" | "indels" => Some(VCF_INDEL << 1),
        "mnp" | "mnps" => Some(VCF_MNP << 1),
        "other" => Some(VCF_OTHER << 1),
        "bnd" => Some(VCF_BND << 1),
        "overlap" => Some(VCF_OVERLAP << 1),
        "ref" => Some(1),
        _ => None,
    }
}

fn extract_gt_query(values: &[Value]) -> Vec<String> {
    let mut out = Vec::new();
    for v in values {
        match v {
            Value::Str(s) => out.push(s.clone()),
            Value::Missing => out.push(".".to_string()),
            _ => {}
        }
    }
    out
}

fn gt_match_eq(gt: &str, query: &[String]) -> bool {
    if query.is_empty() {
        return false;
    }
    for q in query {
        if gt_eq_query(gt, q) {
            return true;
        }
    }
    false
}

fn gt_match_regex(gt: &str, query: &[String], negate: bool) -> bool {
    if query.is_empty() {
        return false;
    }
    for q in query {
        let re = match Regex::new(q) {
            Ok(r) => r,
            Err(_) => continue,
        };
        let m = re.is_match(gt);
        if negate {
            if !m {
                return true;
            }
        } else if m {
            return true;
        }
    }
    false
}

fn gt_eq_query(gt: &str, query: &str) -> bool {
    if query == "." {
        return gt_missing(gt);
    }
    if gt == query {
        return true;
    }
    if query == "aA" || query == "Aa" {
        return gt_is_aA(gt);
    }
    let q = query.to_ascii_lowercase();
    match q.as_str() {
        "rr" => gt_is_rr(gt),
        "ra" | "ar" => gt_is_ra(gt),
        "aa" => gt_is_aa(gt),
        "a" => gt_is_a(gt),
        "r" => gt_is_r(gt),
        "hom" => gt_is_hom(gt),
        "het" => gt_is_het(gt),
        "hap" => gt_is_hap(gt),
        "mis" => gt_missing(gt),
        "ref" => gt_is_ref(gt),
        "alt" => gt_is_alt(gt),
        _ => false,
    }
}

fn gt_missing(gt: &str) -> bool {
    if gt == "." || gt == "./." || gt == ".|." {
        return true;
    }
    gt.split(|c| c == '/' || c == '|')
        .any(|p| p == "." || p.is_empty())
}

fn gt_ploidy(gt: &str) -> usize {
    gt.split(|c| c == '/' || c == '|').count()
}

fn gt_alleles(gt: &str) -> Vec<i32> {
    gt.split(|c| c == '/' || c == '|')
        .filter_map(|p| p.parse::<i32>().ok())
        .collect()
}

fn gt_is_ref(gt: &str) -> bool {
    if gt_missing(gt) {
        return false;
    }
    let alleles = gt_alleles(gt);
    !alleles.is_empty() && alleles.iter().all(|a| *a == 0)
}

fn gt_is_alt(gt: &str) -> bool {
    if gt_missing(gt) {
        return false;
    }
    !gt_is_ref(gt)
}

fn gt_is_rr(gt: &str) -> bool {
    if gt_ploidy(gt) < 2 {
        return false;
    }
    gt_is_ref(gt)
}

fn gt_is_ra(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) < 2 {
        return false;
    }
    let alleles = gt_alleles(gt);
    let has_ref = alleles.iter().any(|a| *a == 0);
    let has_alt = alleles.iter().any(|a| *a > 0);
    has_ref && has_alt
}

fn gt_is_aa(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) < 2 {
        return false;
    }
    let alleles = gt_alleles(gt);
    let all_alt = alleles.iter().all(|a| *a > 0);
    let all_same = alleles.windows(2).all(|w| w[0] == w[1]);
    all_alt && all_same
}

fn gt_is_aA(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) < 2 {
        return false;
    }
    let alleles = gt_alleles(gt);
    let all_alt = alleles.iter().all(|a| *a > 0);
    let all_same = alleles.windows(2).all(|w| w[0] == w[1]);
    all_alt && !all_same
}

fn gt_is_a(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) != 1 {
        return false;
    }
    let alleles = gt_alleles(gt);
    alleles.len() == 1 && alleles[0] > 0
}

fn gt_is_r(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) != 1 {
        return false;
    }
    let alleles = gt_alleles(gt);
    alleles.len() == 1 && alleles[0] == 0
}

fn gt_is_hap(gt: &str) -> bool {
    !gt_missing(gt) && gt_ploidy(gt) == 1
}

fn gt_is_hom(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) < 2 {
        return false;
    }
    let alleles = gt_alleles(gt);
    alleles.windows(2).all(|w| w[0] == w[1])
}

fn gt_is_het(gt: &str) -> bool {
    if gt_missing(gt) || gt_ploidy(gt) < 2 {
        return false;
    }
    let alleles = gt_alleles(gt);
    !alleles.windows(2).all(|w| w[0] == w[1])
}

fn compare_scalar(op: BinaryOp, a: &[Value], b: &[Value]) -> EvalResult {
    let pass = compare_vec(op, a, b);
    EvalResult {
        pass_site: pass,
        pass_samples: None,
    }
}

fn compare_samples_scalar(op: BinaryOp, samples: &SampleValues, other: &[Value]) -> EvalResult {
    let mut pass_samples = vec![false; samples.values.len()];
    let mut pass_site = false;
    for i in 0..samples.values.len() {
        if !samples.mask[i] {
            continue;
        }
        if compare_vec(op, &samples.values[i], other) {
            pass_samples[i] = true;
            pass_site = true;
        }
    }
    EvalResult {
        pass_site,
        pass_samples: Some(pass_samples),
    }
}

fn compare_samples(op: BinaryOp, a: &SampleValues, b: &SampleValues) -> EvalResult {
    let len = a.values.len().min(b.values.len());
    let mut pass_samples = vec![false; len];
    let mut pass_site = false;
    for i in 0..len {
        if !a.mask[i] || !b.mask[i] {
            continue;
        }
        if compare_vec(op, &a.values[i], &b.values[i]) {
            pass_samples[i] = true;
            pass_site = true;
        }
    }
    EvalResult {
        pass_site,
        pass_samples: Some(pass_samples),
    }
}

fn compare_vec(op: BinaryOp, a: &[Value], b: &[Value]) -> bool {
    let miss_one = match op {
        BinaryOp::Eq => false,
        BinaryOp::Ne => true,
        _ => false,
    };
    let miss_two = matches!(op, BinaryOp::Eq);
    if a.is_empty() && b.is_empty() {
        return miss_two;
    }
    if a.is_empty() || b.is_empty() {
        return miss_one;
    }
    for av in a {
        for bv in b {
            if compare_value(op, av, bv, miss_one, miss_two) {
                return true;
            }
        }
    }
    false
}

fn compare_value(op: BinaryOp, a: &Value, b: &Value, miss_one: bool, miss_two: bool) -> bool {
    let a_dot = matches!(a, Value::Str(s) if s == ".");
    let b_dot = matches!(b, Value::Str(s) if s == ".");
    let a_miss = matches!(a, Value::Missing);
    let b_miss = matches!(b, Value::Missing);
    if (a_miss && b_dot) || (b_miss && a_dot) {
        return miss_two;
    }
    if a_miss && b_miss {
        return miss_two;
    }
    if a_miss || b_miss {
        return miss_one;
    }
    match op {
        BinaryOp::Eq => value_eq(a, b),
        BinaryOp::Ne => !value_eq(a, b),
        BinaryOp::Lt => value_cmp(a, b, |x, y| x < y),
        BinaryOp::Le => value_cmp(a, b, |x, y| x <= y),
        BinaryOp::Gt => value_cmp(a, b, |x, y| x > y),
        BinaryOp::Ge => value_cmp(a, b, |x, y| x >= y),
        BinaryOp::Like | BinaryOp::NLike => value_regex(a, b, op == BinaryOp::NLike),
        _ => false,
    }
}

fn value_eq(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Int(x), Value::Int(y)) => x == y,
        (Value::Float(x), Value::Float(y)) => x == y,
        (Value::Int(x), Value::Float(y)) => (*x as f64) == *y,
        (Value::Float(x), Value::Int(y)) => *x == (*y as f64),
        (Value::Str(x), Value::Str(y)) => {
            if x == y {
                true
            } else if x.contains(',') || y.contains(',') {
                let xs: Vec<&str> = x.split(',').collect();
                let ys: Vec<&str> = y.split(',').collect();
                xs.iter().any(|a| ys.iter().any(|b| a == b))
            } else {
                false
            }
        }
        _ => false,
    }
}

fn value_cmp<F: Fn(f64, f64) -> bool>(a: &Value, b: &Value, f: F) -> bool {
    match (a, b) {
        (Value::Int(x), Value::Int(y)) => f(*x as f64, *y as f64),
        (Value::Float(x), Value::Float(y)) => f(*x, *y),
        (Value::Int(x), Value::Float(y)) => f(*x as f64, *y),
        (Value::Float(x), Value::Int(y)) => f(*x, *y as f64),
        _ => false,
    }
}

fn value_regex(a: &Value, b: &Value, negate: bool) -> bool {
    let pattern = match b {
        Value::Str(s) => s,
        _ => return false,
    };
    let re = match Regex::new(pattern) {
        Ok(r) => r,
        Err(_) => return false,
    };
    let target = match a {
        Value::Str(s) => s,
        _ => return false,
    };
    let m = re.is_match(target);
    if negate { !m } else { m }
}

fn value_to_bool(v: EvalValue) -> EvalResult {
    match v {
        EvalValue::Scalar(vs) => EvalResult {
            pass_site: values_truthy(&vs.values),
            pass_samples: None,
        },
        EvalValue::Samples(vs) => {
            let mut pass_samples = vec![false; vs.values.len()];
            let mut pass_site = false;
            for i in 0..vs.values.len() {
                if !vs.mask[i] {
                    continue;
                }
                if values_truthy(&vs.values[i]) {
                    pass_samples[i] = true;
                    pass_site = true;
                }
            }
            EvalResult {
                pass_site,
                pass_samples: Some(pass_samples),
            }
        }
        EvalValue::StrList(_) => EvalResult {
            pass_site: false,
            pass_samples: None,
        },
    }
}

fn values_truthy(values: &[Value]) -> bool {
    for v in values {
        match v {
            Value::Missing => continue,
            Value::Int(i) => {
                if *i != 0 {
                    return true;
                }
            }
            Value::Float(f) => {
                if *f != 0.0 {
                    return true;
                }
            }
            Value::Str(s) => {
                if !s.is_empty() && s != "." {
                    return true;
                }
            }
        }
    }
    false
}

fn apply_unary_neg(v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(mut vs) => {
            for val in vs.values.iter_mut() {
                *val = match val {
                    Value::Int(i) => Value::Int(-*i),
                    Value::Float(f) => Value::Float(-*f),
                    _ => Value::Missing,
                };
            }
            EvalValue::Scalar(vs)
        }
        EvalValue::Samples(mut vs) => {
            for sample in vs.values.iter_mut() {
                for val in sample.iter_mut() {
                    *val = match val {
                        Value::Int(i) => Value::Int(-*i),
                        Value::Float(f) => Value::Float(-*f),
                        _ => Value::Missing,
                    };
                }
            }
            EvalValue::Samples(vs)
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_arith(op: BinaryOp, a: EvalValue, b: EvalValue) -> EvalValue {
    match (a, b) {
        (EvalValue::Scalar(av), EvalValue::Scalar(bv)) => EvalValue::Scalar(ValueVec {
            values: arith_values(op, &av.values, &bv.values),
            is_str: false,
            kind: ValueKind::Normal,
        }),
        (EvalValue::Samples(av), EvalValue::Scalar(bv)) => {
            let mut out = av.values.clone();
            for sample in out.iter_mut() {
                let vals = arith_values(op, sample, &bv.values);
                *sample = vals;
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: av.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (EvalValue::Scalar(av), EvalValue::Samples(bv)) => {
            let mut out = bv.values.clone();
            for sample in out.iter_mut() {
                let vals = arith_values(op, &av.values, sample);
                *sample = vals;
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: bv.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (EvalValue::Samples(av), EvalValue::Samples(bv)) => {
            let len = av.values.len().min(bv.values.len());
            let mut out = vec![Vec::new(); len];
            for i in 0..len {
                out[i] = arith_values(op, &av.values[i], &bv.values[i]);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: av.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (a, _) => a,
    }
}

fn arith_values(op: BinaryOp, a: &[Value], b: &[Value]) -> Vec<Value> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    if a.len() == 1 {
        for bv in b {
            out.push(arith_value(op, &a[0], bv));
        }
        return out;
    }
    if b.len() == 1 {
        for av in a {
            out.push(arith_value(op, av, &b[0]));
        }
        return out;
    }
    let len = a.len().min(b.len());
    for i in 0..len {
        out.push(arith_value(op, &a[i], &b[i]));
    }
    out
}

fn arith_value(op: BinaryOp, a: &Value, b: &Value) -> Value {
    let (x, y) = match (a, b) {
        (Value::Int(x), Value::Int(y)) => (*x as f64, *y as f64),
        (Value::Float(x), Value::Float(y)) => (*x, *y),
        (Value::Int(x), Value::Float(y)) => (*x as f64, *y),
        (Value::Float(x), Value::Int(y)) => (*x, *y as f64),
        _ => return Value::Missing,
    };
    let v = match op {
        BinaryOp::Add => x + y,
        BinaryOp::Sub => x - y,
        BinaryOp::Mul => x * y,
        BinaryOp::Div => x / y,
        BinaryOp::Mod => (x as i64 % y as i64) as f64,
        _ => x,
    };
    Value::Float(v)
}
fn eval_func(name: FuncName, args: &[Expr], ctx: &mut EvalContext) -> Result<EvalValue> {
    match name {
        FuncName::NPass => {
            let expr = args.get(0).ok_or_else(|| anyhow!("N_PASS needs 1 arg"))?;
            let res = eval_bool(expr, ctx)?;
            let count = match res.pass_samples {
                Some(s) => s.iter().filter(|v| **v).count() as i64,
                None => {
                    if res.pass_site {
                        1
                    } else {
                        0
                    }
                }
            };
            Ok(EvalValue::Scalar(ValueVec {
                values: vec![Value::Int(count)],
                is_str: false,
                kind: ValueKind::Normal,
            }))
        }
        FuncName::FPass => {
            let expr = args.get(0).ok_or_else(|| anyhow!("F_PASS needs 1 arg"))?;
            let res = eval_bool(expr, ctx)?;
            let frac = match res.pass_samples {
                Some(s) => {
                    if s.is_empty() {
                        0.0
                    } else {
                        s.iter().filter(|v| **v).count() as f64 / s.len() as f64
                    }
                }
                None => {
                    if res.pass_site {
                        1.0
                    } else {
                        0.0
                    }
                }
            };
            Ok(EvalValue::Scalar(ValueVec {
                values: vec![Value::Float(frac)],
                is_str: false,
                kind: ValueKind::Normal,
            }))
        }
        FuncName::Abs => {
            let v = eval_value(args.get(0).ok_or_else(|| anyhow!("ABS needs 1 arg"))?, ctx)?;
            Ok(apply_abs(v))
        }
        FuncName::Count => {
            let v = eval_value(
                args.get(0).ok_or_else(|| anyhow!("COUNT needs 1 arg"))?,
                ctx,
            )?;
            Ok(apply_count(v))
        }
        FuncName::StrLen => {
            let v = eval_value(
                args.get(0).ok_or_else(|| anyhow!("STRLEN needs 1 arg"))?,
                ctx,
            )?;
            Ok(apply_strlen(v))
        }
        FuncName::Phred => {
            let v = eval_value(
                args.get(0).ok_or_else(|| anyhow!("PHRED needs 1 arg"))?,
                ctx,
            )?;
            Ok(apply_phred(v))
        }
        FuncName::Binom => apply_binom(args, ctx),
        FuncName::Max
        | FuncName::Min
        | FuncName::Avg
        | FuncName::Sum
        | FuncName::Median
        | FuncName::Stdev => {
            let v = eval_value(args.get(0).ok_or_else(|| anyhow!("func needs 1 arg"))?, ctx)?;
            Ok(apply_reduce(name, v))
        }
        FuncName::SMplMax
        | FuncName::SMplMin
        | FuncName::SMplAvg
        | FuncName::SMplSum
        | FuncName::SMplMedian
        | FuncName::SMplStdev
        | FuncName::SMplCount => {
            let v = eval_value(args.get(0).ok_or_else(|| anyhow!("func needs 1 arg"))?, ctx)?;
            Ok(apply_sample_reduce(name, v))
        }
    }
}

fn apply_abs(v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(mut vs) => {
            for val in vs.values.iter_mut() {
                *val = match val {
                    Value::Int(i) => Value::Int(i.abs()),
                    Value::Float(f) => Value::Float(f.abs()),
                    _ => Value::Missing,
                };
            }
            EvalValue::Scalar(vs)
        }
        EvalValue::Samples(mut vs) => {
            for sample in vs.values.iter_mut() {
                for val in sample.iter_mut() {
                    *val = match val {
                        Value::Int(i) => Value::Int(i.abs()),
                        Value::Float(f) => Value::Float(f.abs()),
                        _ => Value::Missing,
                    };
                }
            }
            EvalValue::Samples(vs)
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_count(v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(vs) => {
            let count = vs
                .values
                .iter()
                .filter(|v| !matches!(v, Value::Missing))
                .count() as i64;
            EvalValue::Scalar(ValueVec {
                values: vec![Value::Int(count)],
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::Samples(vs) => {
            let mut count = 0i64;
            for (i, sample) in vs.values.iter().enumerate() {
                if !vs.mask.get(i).copied().unwrap_or(true) {
                    continue;
                }
                count += sample
                    .iter()
                    .filter(|v| !matches!(v, Value::Missing))
                    .count() as i64;
            }
            EvalValue::Scalar(ValueVec {
                values: vec![Value::Int(count)],
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_strlen(v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(vs) => {
            let mut out = Vec::new();
            for val in vs.values {
                match val {
                    Value::Str(s) => out.push(Value::Int(s.len() as i64)),
                    _ => out.push(Value::Missing),
                }
            }
            EvalValue::Scalar(ValueVec {
                values: out,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::Samples(vs) => {
            let mut out = Vec::new();
            for sample in vs.values {
                let mut vals = Vec::new();
                for val in sample {
                    match val {
                        Value::Str(s) => vals.push(Value::Int(s.len() as i64)),
                        _ => vals.push(Value::Missing),
                    }
                }
                out.push(vals);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: vs.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_phred(v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(vs) => {
            let mut out = Vec::new();
            for val in vs.values {
                match val {
                    Value::Float(f) => out.push(Value::Float(-10.0 * f.log10())),
                    Value::Int(i) => out.push(Value::Float(-10.0 * (i as f64).log10())),
                    _ => out.push(Value::Missing),
                }
            }
            EvalValue::Scalar(ValueVec {
                values: out,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::Samples(vs) => {
            let mut out = Vec::new();
            for sample in vs.values {
                let mut vals = Vec::new();
                for val in sample {
                    match val {
                        Value::Float(f) => vals.push(Value::Float(-10.0 * f.log10())),
                        Value::Int(i) => vals.push(Value::Float(-10.0 * (i as f64).log10())),
                        _ => vals.push(Value::Missing),
                    }
                }
                out.push(vals);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: vs.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_binom(args: &[Expr], ctx: &mut EvalContext) -> Result<EvalValue> {
    match args.len() {
        1 => {
            let v = eval_value(&args[0], ctx)?;
            match v {
                EvalValue::Samples(vs) => {
                    let mut out = Vec::with_capacity(vs.values.len());
                    for (i, sample_vals) in vs.values.iter().enumerate() {
                        let (a, b) = if let Some((x, y)) = ad_pair_from_gt(ctx, i, sample_vals) {
                            (x, y)
                        } else {
                            first_two_numbers(sample_vals).unwrap_or((0.0, 0.0))
                        };
                        out.push(vec![Value::Float(binom_two_sided(a, b))]);
                    }
                    Ok(EvalValue::Samples(SampleValues {
                        values: out,
                        mask: vs.mask,
                        is_str: false,
                        kind: ValueKind::Normal,
                    }))
                }
                EvalValue::Scalar(vs) => {
                    let (a, b) = first_two_numbers(&vs.values).unwrap_or((0.0, 0.0));
                    Ok(EvalValue::Scalar(ValueVec {
                        values: vec![Value::Float(binom_two_sided(a, b))],
                        is_str: false,
                        kind: ValueKind::Normal,
                    }))
                }
                EvalValue::StrList(s) => Ok(EvalValue::StrList(s)),
            }
        }
        2 => {
            let a = eval_value(&args[0], ctx)?;
            let b = eval_value(&args[1], ctx)?;
            Ok(apply_binom_two_args(a, b))
        }
        _ => Err(anyhow!("BINOM needs 1 or 2 args")),
    }
}

fn apply_binom_two_args(a: EvalValue, b: EvalValue) -> EvalValue {
    match (a, b) {
        (EvalValue::Scalar(av), EvalValue::Scalar(bv)) => {
            let x = first_number(&av.values).unwrap_or(0.0);
            let y = first_number(&bv.values).unwrap_or(0.0);
            EvalValue::Scalar(ValueVec {
                values: vec![Value::Float(binom_two_sided(x, y))],
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (EvalValue::Samples(av), EvalValue::Scalar(bv)) => {
            let y = first_number(&bv.values).unwrap_or(0.0);
            let mut out = Vec::with_capacity(av.values.len());
            for vals in &av.values {
                let x = first_number(vals).unwrap_or(0.0);
                out.push(vec![Value::Float(binom_two_sided(x, y))]);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: av.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (EvalValue::Scalar(av), EvalValue::Samples(bv)) => {
            let x = first_number(&av.values).unwrap_or(0.0);
            let mut out = Vec::with_capacity(bv.values.len());
            for vals in &bv.values {
                let y = first_number(vals).unwrap_or(0.0);
                out.push(vec![Value::Float(binom_two_sided(x, y))]);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: bv.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (EvalValue::Samples(av), EvalValue::Samples(bv)) => {
            let len = av.values.len().min(bv.values.len());
            let mut out = Vec::with_capacity(len);
            for i in 0..len {
                let x = first_number(&av.values[i]).unwrap_or(0.0);
                let y = first_number(&bv.values[i]).unwrap_or(0.0);
                out.push(vec![Value::Float(binom_two_sided(x, y))]);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: av.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        (a, _) => a,
    }
}

fn first_number(values: &[Value]) -> Option<f64> {
    values.iter().find_map(|v| match v {
        Value::Int(i) => Some(*i as f64),
        Value::Float(f) => Some(*f),
        _ => None,
    })
}

fn first_two_numbers(values: &[Value]) -> Option<(f64, f64)> {
    let mut it = values.iter().filter_map(|v| match v {
        Value::Int(i) => Some(*i as f64),
        Value::Float(f) => Some(*f),
        _ => None,
    });
    let a = it.next()?;
    let b = it.next()?;
    Some((a, b))
}

fn ad_pair_from_gt(
    ctx: &mut EvalContext,
    sample_idx: usize,
    values: &[Value],
) -> Option<(f64, f64)> {
    let gt = parse_gt(ctx, sample_idx);
    if gt.len() < 2 {
        return None;
    }
    let a = gt[0];
    let b = gt[1];
    if a < 0 || b < 0 {
        return None;
    }
    let ai = a as usize;
    let bi = b as usize;
    let get = |i: usize| -> Option<f64> {
        match values.get(i)? {
            Value::Int(v) => Some(*v as f64),
            Value::Float(v) => Some(*v),
            _ => None,
        }
    };
    Some((get(ai)?, get(bi)?))
}

fn binom_two_sided(a: f64, b: f64) -> f64 {
    let x = a.max(0.0).round() as usize;
    let y = b.max(0.0).round() as usize;
    let n = x + y;
    if n == 0 {
        return 0.0;
    }
    let k = x.min(y);
    let mut prob = (0.5f64).powi(n as i32);
    let mut cdf = prob;
    for i in 1..=k {
        prob *= (n - i + 1) as f64 / i as f64;
        cdf += prob;
    }
    (2.0 * cdf).min(1.0)
}

fn apply_reduce(name: FuncName, v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Scalar(vs) => {
            let val = reduce_values(name, &vs.values);
            EvalValue::Scalar(ValueVec {
                values: vec![val],
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::Samples(vs) => {
            let mut all = Vec::new();
            for sample in vs.values {
                all.extend(sample);
            }
            let val = reduce_values(name, &all);
            EvalValue::Scalar(ValueVec {
                values: vec![val],
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        EvalValue::StrList(s) => EvalValue::StrList(s),
    }
}

fn apply_sample_reduce(name: FuncName, v: EvalValue) -> EvalValue {
    match v {
        EvalValue::Samples(vs) => {
            let mut out = Vec::with_capacity(vs.values.len());
            for sample in &vs.values {
                let val = match name {
                    FuncName::SMplCount => Value::Int(
                        sample
                            .iter()
                            .filter(|v| !matches!(v, Value::Missing))
                            .count() as i64,
                    ),
                    _ => reduce_values(sample_reduce_to_func(name), sample),
                };
                out.push(vec![val]);
            }
            EvalValue::Samples(SampleValues {
                values: out,
                mask: vs.mask,
                is_str: false,
                kind: ValueKind::Normal,
            })
        }
        _ => v,
    }
}

fn sample_reduce_to_func(name: FuncName) -> FuncName {
    match name {
        FuncName::SMplMax => FuncName::Max,
        FuncName::SMplMin => FuncName::Min,
        FuncName::SMplAvg => FuncName::Avg,
        FuncName::SMplSum => FuncName::Sum,
        FuncName::SMplMedian => FuncName::Median,
        FuncName::SMplStdev => FuncName::Stdev,
        _ => name,
    }
}

fn reduce_values(name: FuncName, values: &[Value]) -> Value {
    let mut nums: Vec<f64> = values
        .iter()
        .filter_map(|v| match v {
            Value::Int(i) => Some(*i as f64),
            Value::Float(f) => Some(*f),
            _ => None,
        })
        .collect();
    if nums.is_empty() {
        return Value::Missing;
    }
    match name {
        FuncName::Max => Value::Float(nums.iter().cloned().fold(f64::NEG_INFINITY, f64::max)),
        FuncName::Min => Value::Float(nums.iter().cloned().fold(f64::INFINITY, f64::min)),
        FuncName::Sum => Value::Float(nums.iter().sum()),
        FuncName::Avg => Value::Float(nums.iter().sum::<f64>() / nums.len() as f64),
        FuncName::Median => {
            nums.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mid = nums.len() / 2;
            let v = if nums.len() % 2 == 0 {
                (nums[mid - 1] + nums[mid]) / 2.0
            } else {
                nums[mid]
            };
            Value::Float(v)
        }
        FuncName::Stdev => {
            let mean = nums.iter().sum::<f64>() / nums.len() as f64;
            let var = nums.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / nums.len() as f64;
            Value::Float(var.sqrt())
        }
        _ => Value::Missing,
    }
}

fn value_from_bool(res: EvalResult) -> EvalValue {
    if let Some(samples) = res.pass_samples {
        let values = samples
            .iter()
            .map(|b| {
                if *b {
                    vec![Value::Int(1)]
                } else {
                    vec![Value::Missing]
                }
            })
            .collect::<Vec<_>>();
        return EvalValue::Samples(SampleValues {
            values,
            mask: vec![true; samples.len()],
            is_str: false,
            kind: ValueKind::Normal,
        });
    }
    let v = if res.pass_site { 1 } else { 0 };
    EvalValue::Scalar(ValueVec {
        values: vec![Value::Int(v)],
        is_str: false,
        kind: ValueKind::Normal,
    })
}
