# gLasso-f2py

Implement the GLASSOFAST via f2py. Algorithm see:
Sustik, M. A., & Calderhead, B. (2012). GLASSOFAST: an efficient GLASSO implementation. UTCS Technical Report TR-12-29 2012.

Run `python -m f2py.py –c –m glassofast glasso.f90 --fcompiler=gnu95 --compiler=unix` first.
