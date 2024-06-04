*** Copyright Notice ***

“BeamBeam3D”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

 

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

 

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.

 

****************************


BeamBeam3D is a 3D parallel particle-in-cell code for modeling strong-strong or strong-weak beam-beam interactions in high energy ring colliders. This code includes a self-consistent calculation of the electromagnetic forces (i.e. beam-beam forces) from two colliding beams (i.e. strong-strong modeling), a soft-Gaussian approximation of the beam-beam forces, a linear transfer map model for beam transport between collision points, a stochastic map to treat radiation damping, quantum excitation, an arbitrary orbit separation model, a single map to account for chromaticity effects, and models of conducting wire, crab cavity, electron lens for beam-beam compensation. It can handle multiple bunches from each beam collision at multiple interaction points (IPs) with arbitrary separation and crossing angle. The parallel implementation is done using a particle-field decomposition method to achieve a good load balance. It has been applied to studies of the beam-beam effects at a number colliders such as RHIC, Tevatron, LHC, PEP-II, and KEK-B.
