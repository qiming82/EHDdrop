# EHDdrop
electrohydrodynamics of a single drop or drop pairs

This is some test code for a moving boundary problem in classical fluid mechanics. 
The code is set up for a problem for the following particular situations:
  1. fluid part: zero Reynolds number flow (i.e. in Stokes regime) for a viscous drop immersed in another viscous fluids (viscosity \mu_1, \mu_2)
  2. electrostatic part: DC fields with a uniform far field BC, E = -\nabla\phi and E_\infty~ E_0
  3. coupling is realized so-called Taylor-Melcher or leaky dielectric model
Method: boundary integral approach (ref. Pozrikidis 1992 Cambridge text book)
The code can typically handle either
  1. a single viscous drop in uniform electric field (parameters being: viscosity ratio, electric permittivity/conductivity ratios, electric Bond # or Taylor #)
  2. a drop pair with the same or different electric properties
