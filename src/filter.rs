// High-performance VCF filtering using pest grammar (filter_expr.pest)
use crate::vcf::parse_vcf_full_line;
use crate::vcf::VcfParsedRecord;
use anyhow::{anyhow, Result};
use pest::Parser;
use pest_derive::Parser;
use std::collections::HashMap;

#[derive(Parser)]
#[grammar = "filter_expr.pest"]
pub struct ExprParser;

#[derive(Debug, Clone)]
pub enum AstNode {
    And(Box<AstNode>, Box<AstNode>),
    Or(Box<AstNode>, Box<AstNode>),
    Not(Box<AstNode>),
    Cmp {
        left: Box<AstNode>,
        op: CmpOp,
        right: Box<AstNode>,
    },
    Field(String, Option<String>),
    ArrayAccess {
        field: (String, Option<String>),
        index: usize,
    },
    Str(String),
    Int(i64),
    Float(f64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CmpOp {
    Eq,
    Ne,
    Lt,
    Gt,
    Le,
    Ge,
}

#[derive(Debug, Clone)]
pub enum Value {
    Missing,
    Int(i64),
    Float(f64),
    Str(String),
}

impl Value {
    fn cmp(&self, op: CmpOp, other: &Value) -> bool {
        match (self, other) {
            (Value::Missing, _) | (_, Value::Missing) => false,
            (Value::Int(a), Value::Int(b)) => match op {
                CmpOp::Eq => a == b,
                CmpOp::Ne => a != b,
                CmpOp::Lt => a < b,
                CmpOp::Gt => a > b,
                CmpOp::Le => a <= b,
                CmpOp::Ge => a >= b,
            },
            (Value::Float(a), Value::Float(b)) => match op {
                CmpOp::Eq => a == b,
                CmpOp::Ne => a != b,
                CmpOp::Lt => a < b,
                CmpOp::Gt => a > b,
                CmpOp::Le => a <= b,
                CmpOp::Ge => a >= b,
            },
            (Value::Int(a), Value::Float(b)) => self.as_float().cmp(op, &Value::Float(*b)),
            (Value::Float(a), Value::Int(b)) => Value::Float(*a).cmp(op, &Value::Int(*b)),
            (Value::Str(a), Value::Str(b)) => match op {
                CmpOp::Eq => a == b,
                CmpOp::Ne => a != b,
                _ => false,
            },
            _ => false,
        }
    }

    fn as_float(&self) -> Value {
        match self {
            Value::Int(i) => Value::Float(*i as f64),
            v => v.clone(),
        }
    }
}

pub fn eval_ast(ast: &AstNode, rec: &VcfParsedRecord) -> bool {
    match ast {
        AstNode::And(a, b) => eval_ast(a, rec) && eval_ast(b, rec),
        AstNode::Or(a, b) => eval_ast(a, rec) || eval_ast(b, rec),
        AstNode::Not(a) => !eval_ast(a, rec),
        AstNode::Cmp { left, op, right } => {
            let lv = eval_value(left, rec);
            let rv = eval_value(right, rec);
            lv.cmp(*op, &rv)
        }
        _ => false,
    }
}

fn eval_value(node: &AstNode, rec: &VcfParsedRecord) -> Value {
    match node {
        AstNode::Str(s) => Value::Str(s.clone()),
        AstNode::Int(i) => Value::Int(*i),
        AstNode::Float(f) => Value::Float(*f),

        AstNode::Field(key, subkey) => {
            let lookup = if let Some(sub) = subkey {
                format!("{}/{}", key, sub)
            } else {
                key.clone()
            };

            match rec.info.get(&lookup) {
                Some(v_str) => {
                    // Try int
                    if let Ok(i) = v_str.parse::<i64>() {
                        Value::Int(i)
                    } else if let Ok(f) = v_str.parse::<f64>() {
                        Value::Float(f)
                    } else {
                        Value::Str(v_str.clone())
                    }
                }
                None => Value::Missing,
            }
        }

        AstNode::ArrayAccess { field, index } => {
            let lookup = if let Some(sub) = &field.1 {
                format!("{}/{}", field.0, sub)
            } else {
                field.0.clone()
            };
            if let Some(raw) = rec.info.get(&lookup) {
                let parts: Vec<&str> = raw.split(',').collect();
                if *index < parts.len() {
                    let entry = parts[*index];
                    if let Ok(i) = entry.parse::<i64>() {
                        return Value::Int(i);
                    }
                    if let Ok(f) = entry.parse::<f64>() {
                        return Value::Float(f);
                    }
                    return Value::Str(entry.to_string());
                }
            }
            Value::Missing
        }

        _ => Value::Missing,
    }
}

pub fn parse_expr(input: &str) -> Result<AstNode> {
    let pairs = ExprParser::parse(Rule::expr, input)?;
    let ast = build_expr(pairs.into_iter().next().unwrap());
    Ok(ast)
}

fn build_expr(pair: pest::iterators::Pair<Rule>) -> AstNode {
    match pair.as_rule() {
        Rule::or_expr => {
            let mut inner = pair.into_inner();
            let first = build_expr(inner.next().unwrap());
            inner.fold(first, |left, p| {
                let right = build_expr(p.into_inner().next().unwrap());
                AstNode::Or(Box::new(left), Box::new(right))
            })
        }
        Rule::and_expr => {
            let mut inner = pair.into_inner();
            let first = build_expr(inner.next().unwrap());
            inner.fold(first, |left, p| {
                let right = build_expr(p.into_inner().next().unwrap());
                AstNode::And(Box::new(left), Box::new(right))
            })
        }
        Rule::not_expr => {
            let mut inner = pair.into_inner();
            let first = inner.next().unwrap();
            if first.as_rule() == Rule::NOT {
                AstNode::Not(Box::new(build_expr(inner.next().unwrap())))
            } else {
                build_expr(first)
            }
        }
        Rule::cmp_expr => {
            let mut inner = pair.into_inner();
            let left = build_expr(inner.next().unwrap());
            if let Some(op_pair) = inner.next() {
                let op = match op_pair.as_rule() {
                    Rule::EQ => CmpOp::Eq,
                    Rule::NEQ => CmpOp::Ne,
                    Rule::LT => CmpOp::Lt,
                    Rule::GT => CmpOp::Gt,
                    Rule::LE => CmpOp::Le,
                    Rule::GE => CmpOp::Ge,
                    _ => unreachable!(),
                };
                let right = build_expr(inner.next().unwrap());
                AstNode::Cmp {
                    left: Box::new(left),
                    op,
                    right: Box::new(right),
                }
            } else {
                left
            }
        }
        Rule::primary => build_expr(pair.into_inner().next().unwrap()),
        Rule::field => {
            let s = pair.as_str();
            let parts: Vec<&str> = s.split('/').collect();
            if parts.len() == 2 {
                AstNode::Field(parts[0].to_string(), Some(parts[1].to_string()))
            } else {
                AstNode::Field(parts[0].to_string(), None)
            }
        }
        Rule::array_access => {
            let mut inner = pair.into_inner();
            let field_pair = inner.next().unwrap();
            let idx_pair = inner.next().unwrap();
            let s = field_pair.as_str();
            let parts: Vec<&str> = s.split('/').collect();
            let index = idx_pair.as_str().parse::<usize>().unwrap_or(0);
            if parts.len() == 2 {
                AstNode::ArrayAccess {
                    field: (parts[0].to_string(), Some(parts[1].to_string())),
                    index,
                }
            } else {
                AstNode::ArrayAccess {
                    field: (parts[0].to_string(), None),
                    index,
                }
            }
        }
        Rule::string => {
            let s = pair.as_str();
            AstNode::Str(s.trim_matches('"').to_string())
        }
        Rule::int => AstNode::Int(pair.as_str().parse::<i64>().unwrap()),
        Rule::float => AstNode::Float(pair.as_str().parse::<f64>().unwrap()),
        Rule::expr => build_expr(pair.into_inner().next().unwrap()),
        _ => unreachable!("Unhandled rule: {:?}", pair.as_rule()),
    }
}
use crate::filter_args::FilterArgs;

pub fn run_filter(args: &FilterArgs) -> Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    let ast = parse_expr(&args.expr)?;
    let mut reader = crate::vcf::VcfReader::open(&args.input)?;
    let mut out = BufWriter::new(File::create(&args.output)?);

    // Write header unchanged
    for h in &reader.header()? {
        writeln!(out, "{}", h)?;
    }

    let soft = args.soft_filter.clone();
    let pass_only = args.pass_only;

    while let Some((line, _offset)) = reader.next_raw_line()? {
        if let Some(mut rec) = parse_vcf_full_line(&line) {
            let pass = eval_ast(&ast, &rec);

            match (pass, &soft, pass_only) {
                (true, _, false) => rec.filter = "PASS".to_string(),
                (true, _, true) => rec.filter = "PASS".to_string(),
                (false, None, _) => continue,
                (false, Some(name), false) => rec.filter = name.clone(),
                (false, Some(_), true) => continue,
            }

            writeln!(out, "{}", rec.to_line())?;
        }
    }

    Ok(())
}
