#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
