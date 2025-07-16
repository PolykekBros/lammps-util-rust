use std::{
    error, fmt,
    io::{self, BufRead, Lines},
    iter::Enumerate,
    str::FromStr,
};

#[derive(Debug)]
pub struct ParseError {
    line_number: usize,
    source: Box<dyn error::Error>,
}

impl ParseError {
    fn new(line_number: usize, source: Box<dyn error::Error>) -> Self {
        Self {
            line_number,
            source,
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

#[derive(Debug)]
pub struct Token {
    value: String,
    line_index: usize,
}

impl Token {
    fn new(value: String, line_index: usize) -> Self {
        Self { value, line_index }
    }
}

impl TryInto<T> for Token
where
    T: FromStr<Err: error::Error + 'static>,
{
    type Error = ParseError;
    fn try_into(self) -> Result<U, Self::Error> {
        self.value
            .parse::<T>()
            .map_err(|e| ParseError::new(self.line_index, e.into()))
    }
}

#[derive(Debug)]
pub struct TokenIterator<B> {
    lines: Enumerate<Lines<B>>,
    tokens: std::vec::IntoIter<String>,
    line_index: usize,
}

impl<B: BufRead> TokenIterator<B> {
    pub fn new(lines: Lines<B>) -> Self {
        Self {
            lines: lines.enumerate(),
            tokens: Default::default(),
            line_index: 0,
        }
    }
}

impl<B: BufRead> Iterator for TokenIterator<B> {
    type Item = Result<Token, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.tokens.next() {
            Some(s) => Some(Ok(Token::new(s, self.line_index))),
            None => match self.lines.next()? {
                (i, Ok(s)) => {
                    self.line_index = 1;
                    self.tokens = s
                        .split_whitespace()
                        .map(str::to_string)
                        .collect::<Vec<_>>()
                        .into_iter();
                    self.next()
                }
                (i, Err(e)) => Some(Err(ParseError::new(i, e.into()))),
            },
        }
    }
}

pub struct Parser<B> {
    iter: TokenIterator<B>,
}

impl<B: BufRead> Parser<B> {
    pub fn new(lines: Lines<B>) -> Self {
        Self {
            iter: TokenIterator::new(lines),
        }
    }

    pub fn next<T>(&mut self) -> Option<Result<T, ParseError>>
    where
        T: FromStr<Err: error::Error + 'static>,
    {
        Some(
            self.iter
                .next()?
                .and_then(|(i, s)| s.parse::<T>().map_err(|e| ParseError::new(i, e.into()))),
        )
    }
}

fn token_iterator(
    lines: impl IntoIterator<Item = io::Result<String>>,
) -> impl Iterator<Item = Result<String, ParseError>> {
    lines.into_iter().enumerate().map(|(i, r)| (i + 1, r))
}
