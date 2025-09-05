//
// Created by jannickb on 10/2/24.
//
#pragma once

#include <iostream>
#include <kamping/communicator.hpp>

#define KADISREDU_RLOG kamping::comm_world().is_root() && kadisredu::Logger() << "[LOG]: "
#define KADISREDU_RLOG_ERROR kamping::comm_world().is_root() && kadisredu::Logger() << "[ERROR]: "
#define KADISREDU_LOG kadisredu::Logger() << "[LOG][" << kamping::comm_world().rank() << "]: "
#define KADISREDU_LOG_ERROR kadisredu::Logger() << "[ERROR][" << kamping::comm_world().rank() << "]: "
#define KADISREDU_SLOG kadisredu::SyncLogger() << "[LOG][" << kamping::comm_world().rank() << "]: "

// This logger class is based on the one by Daniel Seemaier (KaMinPar)
namespace kadisredu{

class Logger {
 public:


  explicit Logger(std::ostream& out = std::cout, std::string_view end = "\n")
      : _out(out), _end(end) {}

  void flush() {
    if (!flushed) {
      _out << _buf.str() << _end << std::flush;
      flushed = true;
    }
  }

  template <typename Arg>
  Logger& operator<<(Arg&& arg) {
    _buf << arg;
    return *this;
  }

  ~Logger() { flush(); }

  operator bool() {
    return false;
  }

 protected:
  std::stringstream _buf;
  std::ostream& _out;
  std::string_view _end;
  bool flushed{false};
};

class SyncLogger {
 public:
  template <typename... Args>
  explicit SyncLogger(MPI_Comm comm = MPI_COMM_WORLD, std::ostream& out = std::cout, std::string_view end = "\n") : _logger(out, end), _out(out), _end(end), _comm(comm) {};

  void flush() {
    _comm.barrier();
    for (std::size_t rank = 0; rank < _comm.size(); ++rank) {

      if (rank == _comm.rank()) {
        _out << "------- PE=" << rank << " -------" << _end;
        _logger.flush();
        _out << "--------------" << _end;
        _out.flush();
      }
      _comm.barrier();
    }
  }

  ~SyncLogger() { flush(); }

  operator bool() {
    return false;
  }

  template <typename Arg>
  SyncLogger& operator<<(Arg&& arg) {
    _logger << arg;
    return *this;
  }

 protected:
  Logger _logger;
  std::ostringstream _buf;
  std::ostream& _out;
  std::string_view _end;
  kamping::Communicator<> _comm;
};

}