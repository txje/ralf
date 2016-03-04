extern crate log;

use log::{LogRecord, LogLevel, LogMetadata};
use std::io;
use std::io::Write;

struct EasyLogger;

/*
pub enum LogLevel {
  Error, // most serious
  Warn,
  Info,
  Debug,
  Trace, // most verbose
}
*/

impl log::Log for EasyLogger {
  fn enabled(&self, metadata: &LogMetadata) -> bool {
    metadata.level() >= LogLevel::Error
  }

  fn log(&self, record: &LogRecord) {
    if self.enabled(record.metadata()) {
      let _ = match record.level() {
        LogLevel::Trace => writeln!(&mut io::stderr(), "{} - {}", record.level(), record.args()),
        LogLevel::Debug => writeln!(&mut io::stderr(), "{} - {}", record.level(), record.args()),
        LogLevel::Info => writeln!(&mut io::stderr(), "{} - {}", record.level(), record.args()),
        LogLevel::Warn => writeln!(&mut io::stderr(), "{} - {}", record.level(), record.args()),
        LogLevel::Error => writeln!(&mut io::stderr(), "{} - {}", record.level(), record.args()),
      };
    }
  }
}

pub fn init() -> Result<(), log::SetLoggerError> {
  log::set_logger(|max_log_level| {
    max_log_level.set(log::LogLevelFilter::Trace);
    Box::new(EasyLogger)
  })
}
