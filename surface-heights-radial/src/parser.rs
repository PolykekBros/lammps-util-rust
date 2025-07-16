use std::{
    error, fmt,
    io::{self},
    iter::{self},
    str,
};

#[derive(Debug)]
pub struct ParseError {
    line_number: usize,
    source: Box<dyn error::Error>,
}

impl ParseError {
    fn new(line_number: usize, source: impl error::Error + 'static) -> Self {
        Self {
            line_number,
            source: source.into(),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line_number, self.source)
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        Some(&*self.source)
    }
}

type ParseResult<T> = Result<T, ParseError>;

#[derive(Debug)]
pub struct Token {
    value: String,
    line_number: usize,
}

impl Token {
    fn new(value: String, line_number: usize) -> Self {
        Self { value, line_number }
    }
}

impl<'a> From<(usize, &'a str)> for Token {
    fn from(value: (usize, &'a str)) -> Self {
        Token::new(value.1.into(), value.0)
    }
}

fn token_iterator(
    lines: impl IntoIterator<Item = io::Result<String>>,
) -> impl Iterator<Item = ParseResult<Token>> {
    iter::zip(1.., lines).flat_map(|(line_number, str)| match str {
        Ok(str) => str
            .split_whitespace()
            .map(|s| Ok((line_number, s).into()))
            .collect::<Vec<_>>(),
        Err(err) => vec![Err(ParseError::new(line_number, err))],
    })
}

pub struct Parser {
    iter: Box<dyn Iterator<Item = ParseResult<Token>>>,
}

impl Parser {
    pub fn new(lines: impl IntoIterator<Item = io::Result<String>> + 'static) -> Self {
        Self {
            iter: Box::new(token_iterator(lines)),
        }
    }

    pub fn next_parse<T>(&mut self) -> Option<ParseResult<T>>
    where
        T: str::FromStr<Err: error::Error + 'static>,
    {
        self.iter.next().map(|token| {
            token.and_then(|token| {
                token
                    .value
                    .parse::<T>()
                    .map_err(|err| ParseError::new(token.line_number, err))
            })
        })
    }
}

impl Iterator for Parser {
    type Item = ParseResult<String>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|token| token.map(|token| token.value))
    }
}
