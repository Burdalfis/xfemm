-- One-point closure check for the persistent AGE cycle benchmark.
showconsole()
open("base_motor.fem")
DEG = 3.141592653589793 / 180
angle = 360
theta_line = 6 * angle - 30
line_a = sin(theta_line * DEG)
line_b = sin((theta_line - 120) * DEG)
line_c = sin((theta_line + 120) * DEG)
mi_modifycircprop("A", 1, (line_a - line_b) / 3)
mi_modifycircprop("B", 1, (line_b - line_c) / 3)
mi_modifycircprop("C", 1, (line_c - line_a) / 3)
mi_modifyboundprop("SlidingBand", 11, angle)
mi_saveas("age_360.fem")
mi_analyze(1)
mi_loadsolution()
torque = mo_gapintegral("SlidingBand", 0)
ia, va, lambda_a = mo_getcircuitproperties("A")
ib, vb, lambda_b = mo_getcircuitproperties("B")
ic, vc, lambda_c = mo_getcircuitproperties("C")
gap_a_0 = mo_getgapa("SlidingBand", 0)
h1_acc, h1_acs, h1_brc, h1_brs, h1_btc, h1_bts =
    mo_getgapharmonics("SlidingBand", 1)
h6_acc, h6_acs, h6_brc, h6_brs, h6_btc, h6_bts =
    mo_getgapharmonics("SlidingBand", 6)
print("XFEMM_AGE_CLOSURE_RESULT," .. angle .. "," .. torque .. "," ..
      lambda_a .. "," .. lambda_b .. "," .. lambda_c .. "," .. gap_a_0 ..
      "," .. h1_brc .. "," .. h1_brs .. "," .. h1_btc .. "," .. h1_bts ..
      "," .. h6_brc .. "," .. h6_brs .. "," .. h6_btc .. "," .. h6_bts)
mo_close()
quit()
