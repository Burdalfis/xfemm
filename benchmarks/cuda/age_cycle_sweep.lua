-- Persistent sliding-air-gap topology/physics sweep.
-- Run in a directory containing base_motor.fem.  The model is loaded once;
-- mesh reuse and warm starting are controlled by the normal XFEMM sweep
-- environment variables.

showconsole()
open("base_motor.fem")

DEG = 3.141592653589793 / 180
Npoles = 12
current_peak = 1.0

function set_balanced_delta_iq(Iq_line_peak, theta_mech)
    theta_line = (Npoles/2) * theta_mech - 30
    line_a = Iq_line_peak * sin(theta_line * DEG)
    line_b = Iq_line_peak * sin((theta_line - 120) * DEG)
    line_c = Iq_line_peak * sin((theta_line + 120) * DEG)
    wind_a = (line_a - line_b) / 3
    wind_b = (line_b - line_c) / 3
    wind_c = (line_c - line_a) / 3
    mi_modifycircprop("A", 1, wind_a)
    mi_modifycircprop("B", 1, wind_b)
    mi_modifycircprop("C", 1, wind_c)
end

print("XFEMM_AGE_CYCLE_HEADER,angle,current,torque,lambda_a,lambda_b,lambda_c," ..
      "gap_a_0,h1_brc,h1_brs,h1_btc,h1_bts,h6_brc,h6_brs,h6_btc,h6_bts")

for angle = 0, 360, 5 do
    set_balanced_delta_iq(current_peak, angle)
    mi_modifyboundprop("SlidingBand", 11, angle)
    mi_saveas("age_" .. angle .. ".fem")
    print("XFEMM_AGE_CYCLE_BEGIN," .. angle)
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

    print("XFEMM_AGE_CYCLE_RESULT," .. angle .. "," .. current_peak .. "," ..
          torque .. "," .. lambda_a .. "," .. lambda_b .. "," .. lambda_c ..
          "," .. gap_a_0 .. "," .. h1_brc .. "," .. h1_brs .. "," ..
          h1_btc .. "," .. h1_bts .. "," .. h6_brc .. "," .. h6_brs ..
          "," .. h6_btc .. "," .. h6_bts)
    mo_close()
end

quit()
