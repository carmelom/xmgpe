# A working plan

phisical input:
- mass (atomic units)
- scatt length (Borh)
- a transverse trap_freq if ndim < 3 to renormalize the scatt_len
- one to three dimensions with
    - trap_freq
    - geometry
    - domain:
        - if tuple, use it
        - if single number is right boundary, and compute left boundary from geometry
        - if startswith("x"), is a multiple of thomasfermi radius
    - lattice
        - if single number, use it
        - if startswith("x") compute such that the spacing is multiple of healing length

computational parameters:
- 
