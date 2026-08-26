# Third-party notices

The UVLM backend is derived from
[VortexLattice.jl](https://github.com/byuflowlab/VortexLattice.jl), originally
authored by Taylor McDonnell, Andrew Ning, and contributors. VortexLattice.jl is
distributed under the MIT License. The corresponding license text is retained
in `LICENSE`.

The imported backend contains later wing--propeller research modifications from
[`Wing_Propeller_UVLM`](https://github.com/jptpsantos/Wing_Propeller_UVLM).

The near-field force implementation in `src/backend/nearfield.jl` is a
Julia adaptation of the force equations and force-to-node mapping in
[`ImperialCollegeLondon/UVLM`](https://github.com/ImperialCollegeLondon/UVLM),
revision `d8af34a22baf1cddd38f1e362274c407637aab1c`. That project is distributed
under the following BSD 3-Clause License:

Copyright (c) 2018, Imperial College London
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
