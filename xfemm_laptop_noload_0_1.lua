-- ============================================================
--  1404 outrunner - geometry v2
--  9 slots / 12 poles
--  Dimensions in mm
--
--  IMPORTANT:
--    FEMM Lua trig functions use RADIANS.
--    FEMM mi_addarc() angle uses DEGREES.
--  Internally, geometric angles below are kept in degrees,
--  and converted to radians only when calling trig functions.
-- ============================================================

showconsole()
newdocument(0)

mi_probdef(0, "millimeters", "planar", 1e-8, 4.05, 30)

DEG = 3.141592653589793 / 180


-- ============================================================
-- USER DIMENSIONS
-- ============================================================

Nslots = 9
Npoles = 12

iron_corner_r = 0.2
mag_corner_r  = 0.1

-- Stator
R_bore       = 6.1 / 2
R_slotbottom = 7.9 / 2
R_stator     = 14.05 / 2

tooth_body_w = 1.2
tooth_tip_w  = 3.75
tip_radial   = 0.9

-- Rotor / magnets
R_mag_inner  = 14.35 / 2
mag_thick    = 1
mag_width    = 3.3

-- Approximate back iron
back_thick   = 1

-- Outer simulation boundary
R_boundary   = 25

-- ============================================================
-- DERIVED DIMENSIONS
-- ============================================================

slot_pitch = 360 / Nslots     -- 40 deg
pole_pitch = 360 / Npoles     -- 30 deg

half_body = tooth_body_w / 2
half_tip  = tooth_tip_w  / 2

R_tip_inner = R_stator - tip_radial
R_mag_outer = R_mag_inner + mag_thick

-- Back iron starts directly behind the constant-thickness magnets.
R_back_inner = R_mag_outer
R_back_outer = R_back_inner + back_thick

-- asin() returns RADIANS, so convert these back to degrees.
a_root = asin(half_body / R_slotbottom) / DEG
a_body = asin(half_body / R_tip_inner)  / DEG
a_tip  = asin(half_tip  / R_stator)     / DEG

print("Air gap                  = ", R_mag_inner - R_stator, " mm")
print("Stator yoke thickness    = ", R_slotbottom - R_bore, " mm")
print("Slot pitch               = ", slot_pitch, " deg")
print("Pole pitch               = ", pole_pitch, " deg")
print("Tooth-tip total angle    = ", 2*a_tip, " deg")
print("Slot opening angle       = ", slot_pitch - 2*a_tip, " deg")


-- ============================================================
-- HELPERS
-- ============================================================

-- Input angle in DEGREES.
function polar(r, a)
    ar = a * DEG
    return r*cos(ar), r*sin(ar)
end

-- Rotate point by angle in DEGREES.
function rotatexy(x, y, a)
    ar = a * DEG
    xr = x*cos(ar) - y*sin(ar)
    yr = x*sin(ar) + y*cos(ar)
    return xr, yr
end

function addcircle(r)
    mi_addnode( r, 0)
    mi_addnode(-r, 0)

    -- mi_addarc angle is in DEGREES
    mi_addarc( r, 0, -r, 0, 180, 3)
    mi_addarc(-r, 0,  r, 0, 180, 3)
end


-- ============================================================
-- STATOR BORE
-- ============================================================

addcircle(R_bore)


tooth_shoe_profile = "straight_inner"

if tooth_shoe_profile == "straight_inner" then

-- ============================================================
-- STATOR OUTER PROFILE - V3
--
-- Tooth shoe is modeled as a true T-shape:
--   * 0.90 mm radial thickness at tooth centreline
--   * straight, colinear inner shoe faces
--   * 0.30 mm fillets on all stator corners
--   * 0.97 mm FINAL minimum slot opening between the
--     opposing inner shoe-end fillet arcs
--   * circular R_slotbottom arcs between teeth
-- ============================================================

slot_opening = 1

half_body = tooth_body_w / 2

-- This is the inner face of the tooth shoe.
-- Since the stator OD is R_stator, this makes the shoe exactly
-- tip_radial thick at the tooth centreline.
X_tip_inner = R_stator - tip_radial

-- ------------------------------------------------------------
-- Tooth-root / slot-bottom fillet
--
-- Upper root is constructed explicitly. Lower root is its mirror.
-- The fillet is tangent to:
--   y = +half_body
-- and
--   radius R_slotbottom
-- ------------------------------------------------------------

root_fillet_cy = half_body + iron_corner_r
root_fillet_cr = R_slotbottom + iron_corner_r

root_fillet_cx =
    sqrt(
        root_fillet_cr*root_fillet_cr
        - root_fillet_cy*root_fillet_cy
    )

-- Tangency on upper straight tooth-body side
RLx = root_fillet_cx
RLy = half_body

-- Tangency on slot-bottom circle
RSx = root_fillet_cx * R_slotbottom / root_fillet_cr
RSy = root_fillet_cy * R_slotbottom / root_fillet_cr

root_slot_angle =
    atan2(RSy, RSx) / DEG

-- Arc angle around root fillet centre:
-- vector to RS -> vector to RL
rv1x = RSx - root_fillet_cx
rv1y = RSy - root_fillet_cy

rv2x = RLx - root_fillet_cx
rv2y = RLy - root_fillet_cy

rcang =
    (rv1x*rv2x + rv1y*rv2y)
    / (iron_corner_r*iron_corner_r)

if rcang > 1 then rcang = 1 end
if rcang < -1 then rcang = -1 end

root_fillet_angle = acos(rcang) / DEG

-- The actual circular yoke OD remaining between adjacent
-- root fillets.
slot_bottom_angle =
    slot_pitch - 2*root_slot_angle


-- ------------------------------------------------------------
-- Body-to-shoe concave fillet
--
-- Sharp theoretical corner:
--   (X_tip_inner, half_body)
--
-- The 0.30 mm fillet is tangent to the straight tooth body
-- and the straight inner face of the shoe.
-- ------------------------------------------------------------

BTx = X_tip_inner - iron_corner_r
BTy = half_body

BIx = X_tip_inner
BIy = half_body + iron_corner_r

-- This corner is exactly 90 degrees.
body_tip_fillet_angle = 90


-- ------------------------------------------------------------
-- Shoe circumferential extent
--
-- The inner shoe-end fillet is the feature that gets closest
-- to the neighbouring tooth. Solve its Y position so the
-- finished minimum gap between the two opposing r=0.30 mm
-- fillet arcs is exactly slot_opening.
--
-- Slot centreline is half a slot pitch from tooth centre.
-- ------------------------------------------------------------

slot_half_angle = (slot_pitch / 2) * DEG

Y_tip_end =
    iron_corner_r
    + (
        sin(slot_half_angle)
        * (X_tip_inner + iron_corner_r)
        - iron_corner_r
        - slot_opening/2
      )
      / cos(slot_half_angle)

-- Inner shoe-end fillet:
-- centre = (X_tip_inner+r, Y_tip_end-r)
--
-- Tangency to straight inner shoe face:
IEx = X_tip_inner
IEy = Y_tip_end - iron_corner_r

-- Tangency to straight end face:
ETx = X_tip_inner + iron_corner_r
ETy = Y_tip_end


-- ------------------------------------------------------------
-- Outer shoe-end fillet
--
-- The tooth end face is straight and parallel to the local
-- radial X axis. Its outer 0.30 mm fillet is tangent to the
-- stator OD circle.
-- ------------------------------------------------------------

outer_end_fillet_cy = Y_tip_end - iron_corner_r
outer_end_fillet_cr = R_stator - iron_corner_r

outer_end_fillet_cx =
    sqrt(
        outer_end_fillet_cr*outer_end_fillet_cr
        - outer_end_fillet_cy*outer_end_fillet_cy
    )

-- Tangency to straight end face
OTx = outer_end_fillet_cx
OTy = Y_tip_end

-- Tangency to stator OD
OOx =
    outer_end_fillet_cx
    * R_stator
    / outer_end_fillet_cr

OOy =
    outer_end_fillet_cy
    * R_stator
    / outer_end_fillet_cr

outer_face_half_angle =
    atan2(OOy, OOx) / DEG

-- Around the outer fillet centre, the end-face tangent is at
-- 90 degrees and the OD tangent is at outer_face_half_angle.
outer_end_fillet_angle =
    90 - outer_face_half_angle


-- ------------------------------------------------------------
-- Diagnostics / dimensional checks
-- ------------------------------------------------------------

modeled_tip_max_width = 2 * Y_tip_end
modeled_tip_center_thickness = R_stator - X_tip_inner

print("Final slot opening           = ", slot_opening, " mm")
print("Tooth shoe centre thickness  = ", modeled_tip_center_thickness, " mm")
print("Modeled max tooth-tip width  = ", modeled_tip_max_width, " mm")
print("Measured rough tooth width   = ", tooth_tip_w, " mm")
print("Root fillet arc angle        = ", root_fillet_angle, " deg")
print("Slot-bottom arc angle        = ", slot_bottom_angle, " deg")
print("Outer shoe fillet arc angle  = ", outer_end_fillet_angle, " deg")


-- ============================================================
-- DRAW ALL 9 TEETH
-- ============================================================

for i = 0, Nslots-1 do

    t = i * slot_pitch

    -- ---------------- LOWER HALF ----------------

    -- Root fillet tangencies
    RSlx, RSly = rotatexy(RSx, -RSy, t)
    RLlx, RLly = rotatexy(RLx, -RLy, t)

    -- Tooth body / inner shoe fillet
    BTlx, BTly = rotatexy(BTx, -BTy, t)
    BIlx, BIly = rotatexy(BIx, -BIy, t)

    -- Inner shoe-end fillet
    IElx, IEly = rotatexy(IEx, -IEy, t)
    ETlx, ETly = rotatexy(ETx, -ETy, t)

    -- Outer shoe-end fillet
    OTlx, OTly = rotatexy(OTx, -OTy, t)
    OOlx, OOly = rotatexy(OOx, -OOy, t)


    -- ---------------- UPPER HALF ----------------

    OOUx, OOUy = rotatexy(OOx, OOy, t)
    OTUx, OTUy = rotatexy(OTx, OTy, t)

    ETUx, ETUy = rotatexy(ETx, ETy, t)
    IEUx, IEUy = rotatexy(IEx, IEy, t)

    BIUx, BIUy = rotatexy(BIx, BIy, t)
    BTUx, BTUy = rotatexy(BTx, BTy, t)

    RLUx, RLUy = rotatexy(RLx, RLy, t)
    RSUx, RSUy = rotatexy(RSx, RSy, t)


    -- Add nodes
    mi_addnode(RSlx, RSly)
    mi_addnode(RLlx, RLly)

    mi_addnode(BTlx, BTly)
    mi_addnode(BIlx, BIly)

    mi_addnode(IElx, IEly)
    mi_addnode(ETlx, ETly)

    mi_addnode(OTlx, OTly)
    mi_addnode(OOlx, OOly)

    mi_addnode(OOUx, OOUy)
    mi_addnode(OTUx, OTUy)

    mi_addnode(ETUx, ETUy)
    mi_addnode(IEUx, IEUy)

    mi_addnode(BIUx, BIUy)
    mi_addnode(BTUx, BTUy)

    mi_addnode(RLUx, RLUy)
    mi_addnode(RSUx, RSUy)


    -- --------------------------------------------------------
    -- LOWER boundary from slot bottom outward
    -- --------------------------------------------------------

    -- Lower root fillet
    mi_addarc(
        RLlx, RLly,
        RSlx, RSly,
        root_fillet_angle,
        1
    )

    -- Lower straight tooth body
    mi_addsegment(
        RLlx, RLly,
        BTlx, BTly
    )

    -- Lower body-to-shoe fillet
    mi_addarc(
        BIlx, BIly,
        BTlx, BTly,
        body_tip_fillet_angle,
        1
    )

    -- Lower straight INNER SHOE FACE
    mi_addsegment(
        BIlx, BIly,
        IElx, IEly
    )

    -- Lower inner shoe-end fillet
    mi_addarc(
        IElx, IEly,
        ETlx, ETly,
        90,
        1
    )

    -- Lower straight shoe end face
    mi_addsegment(
        ETlx, ETly,
        OTlx, OTly
    )

    -- Lower outer shoe-end fillet
    mi_addarc(
        OTlx, OTly,
        OOlx, OOly,
        outer_end_fillet_angle,
        1
    )


    -- --------------------------------------------------------
    -- Curved OD face of tooth shoe
    -- --------------------------------------------------------

    mi_addarc(
        OOlx, OOly,
        OOUx, OOUy,
        2*outer_face_half_angle,
        1
    )


    -- --------------------------------------------------------
    -- UPPER boundary back inward
    -- --------------------------------------------------------

    -- Upper outer shoe-end fillet
    mi_addarc(
        OOUx, OOUy,
        OTUx, OTUy,
        outer_end_fillet_angle,
        1
    )

    -- Upper straight shoe end face
    mi_addsegment(
        OTUx, OTUy,
        ETUx, ETUy
    )

    -- Upper inner shoe-end fillet
    mi_addarc(
        ETUx, ETUy,
        IEUx, IEUy,
        90,
        1
    )

    -- Upper straight INNER SHOE FACE
    mi_addsegment(
        IEUx, IEUy,
        BIUx, BIUy
    )

    -- Upper body-to-shoe fillet
    mi_addarc(
        BTUx, BTUy,
        BIUx, BIUy,
        body_tip_fillet_angle,
        1
    )

    -- Upper straight tooth body
    mi_addsegment(
        BTUx, BTUy,
        RLUx, RLUy
    )

    -- Upper root fillet
    mi_addarc(
        RSUx, RSUy,
        RLUx, RLUy,
        root_fillet_angle,
        1
    )


    -- --------------------------------------------------------
    -- EXPLICIT circular slot-bottom / stator-yoke OD arc
    -- from this tooth's upper root fillet to the NEXT tooth's
    -- lower root fillet.
    -- --------------------------------------------------------

    next_t = (i + 1) * slot_pitch

    NextRSlx, NextRSly =
        rotatexy(RSx, -RSy, next_t)

    mi_addnode(NextRSlx, NextRSly)

    mi_addarc(
        RSUx, RSUy,
        NextRSlx, NextRSly,
        slot_bottom_angle,
        1
    )
end



elseif tooth_shoe_profile == "concentric_inner" then

-- ============================================================
-- CONCENTRIC-INNER TOOTH SHOE
--
-- Inner and outer tooth-shoe faces are concentric circular arcs.
-- The shoe therefore has constant radial thickness before corner fillets.
--
-- The final slot opening is solved between the opposing INNER
-- shoe-end fillet arcs, analogous to the straight-inner model.
-- ============================================================

slot_opening = 1
half_body = tooth_body_w / 2

-- ------------------------------------------------------------
-- Tooth-root / slot-bottom fillet
-- ------------------------------------------------------------

root_fillet_cy = half_body + iron_corner_r
root_fillet_cr = R_slotbottom + iron_corner_r

root_fillet_cx =
    sqrt(
        root_fillet_cr*root_fillet_cr
        - root_fillet_cy*root_fillet_cy
    )

RLx = root_fillet_cx
RLy = half_body

RSx = root_fillet_cx * R_slotbottom / root_fillet_cr
RSy = root_fillet_cy * R_slotbottom / root_fillet_cr

root_slot_angle = atan2(RSy, RSx) / DEG

rv1x = RSx - root_fillet_cx
rv1y = RSy - root_fillet_cy
rv2x = RLx - root_fillet_cx
rv2y = RLy - root_fillet_cy

rcang =
    (rv1x*rv2x + rv1y*rv2y)
    / (iron_corner_r*iron_corner_r)

if rcang > 1 then rcang = 1 end
if rcang < -1 then rcang = -1 end

root_fillet_angle = acos(rcang) / DEG
slot_bottom_angle = slot_pitch - 2*root_slot_angle

-- ------------------------------------------------------------
-- Body-to-inner-arc concave fillet
--
-- Fillet centre is tangent to:
--   y = half_body
-- and the concentric inner shoe circle R_tip_inner.
-- ------------------------------------------------------------

body_fillet_center_r = R_tip_inner - iron_corner_r
body_fillet_cy = half_body + iron_corner_r

if body_fillet_center_r <= body_fillet_cy then
    error("Concentric shoe: body fillet geometry is impossible")
end

body_fillet_cx =
    sqrt(
        body_fillet_center_r*body_fillet_center_r
        - body_fillet_cy*body_fillet_cy
    )

BTx = body_fillet_cx
BTy = half_body

-- Tangency to inner shoe circle.
BIx = body_fillet_cx * R_tip_inner / body_fillet_center_r
BIy = body_fillet_cy * R_tip_inner / body_fillet_center_r

body_fillet_center_angle =
    atan2(body_fillet_cy, body_fillet_cx) / DEG

body_inner_tangent_angle =
    atan2(BIy, BIx) / DEG

body_tip_fillet_angle =
    90 + body_fillet_center_angle

-- ------------------------------------------------------------
-- Inner shoe-end fillet and final slot opening
--
-- The inner fillet is tangent to:
--   * inner concentric circle R_tip_inner
--   * a radial tooth-end face
--
-- Its centre lies at radius R_tip_inner + r.
-- Solve its centre angle so the gap between the two opposing
-- equal-radius fillet circles is exactly slot_opening.
-- ------------------------------------------------------------

inner_end_center_r = R_tip_inner + iron_corner_r

gap_arg =
    (slot_opening + 2*iron_corner_r)
    / (2*inner_end_center_r)

if gap_arg >= 1 then
    error("Concentric shoe: requested slot opening is too large")
end

inner_end_center_angle =
    slot_pitch/2
    - asin(gap_arg) / DEG

inner_end_delta =
    asin(
        iron_corner_r
        / inner_end_center_r
    ) / DEG

-- Radial end-face angle.
end_face_half_angle =
    inner_end_center_angle
    + inner_end_delta

-- Tangency on inner concentric shoe face.
IEx, IEy =
    polar(
        R_tip_inner,
        inner_end_center_angle
    )

-- Tangency on radial end face.
inner_end_tangent_r =
    inner_end_center_r
    * cos(inner_end_delta*DEG)

ETx, ETy =
    polar(
        inner_end_tangent_r,
        end_face_half_angle
    )

inner_end_fillet_angle =
    90 - inner_end_delta

-- ------------------------------------------------------------
-- Outer shoe-end fillet
--
-- Tangent to:
--   * stator OD circle R_stator
--   * same radial tooth-end face
-- ------------------------------------------------------------

outer_end_center_r =
    R_stator - iron_corner_r

outer_end_delta =
    asin(
        iron_corner_r
        / outer_end_center_r
    ) / DEG

outer_end_center_angle =
    end_face_half_angle
    - outer_end_delta

-- Tangency on OD circle.
OOx, OOy =
    polar(
        R_stator,
        outer_end_center_angle
    )

-- Tangency on radial end face.
outer_end_tangent_r =
    outer_end_center_r
    * cos(outer_end_delta*DEG)

OTx, OTy =
    polar(
        outer_end_tangent_r,
        end_face_half_angle
    )

outer_end_fillet_angle =
    90 + outer_end_delta

outer_face_half_angle =
    outer_end_center_angle

inner_shoe_arc_angle =
    inner_end_center_angle
    - body_inner_tangent_angle

modeled_tip_max_width =
    2 * OTy

print("Tooth shoe profile          = concentric_inner")
print("Final slot opening          = ", slot_opening, " mm")
print("Inner shoe radius           = ", R_tip_inner, " mm")
print("Outer shoe radius           = ", R_stator, " mm")
print("Shoe radial thickness       = ", R_stator - R_tip_inner, " mm")
print("Modeled max tooth-tip width = ", modeled_tip_max_width, " mm")
print("Measured rough tooth width  = ", tooth_tip_w, " mm")
print("Radial end-face half-angle  = ", end_face_half_angle, " deg")

-- ============================================================
-- DRAW ALL TEETH
-- ============================================================

for i = 0, Nslots-1 do

    t = i * slot_pitch

    -- LOWER HALF
    RSlx, RSly = rotatexy(RSx, -RSy, t)
    RLlx, RLly = rotatexy(RLx, -RLy, t)

    BTlx, BTly = rotatexy(BTx, -BTy, t)
    BIlx, BIly = rotatexy(BIx, -BIy, t)

    IElx, IEly = rotatexy(IEx, -IEy, t)
    ETlx, ETly = rotatexy(ETx, -ETy, t)

    OTlx, OTly = rotatexy(OTx, -OTy, t)
    OOlx, OOly = rotatexy(OOx, -OOy, t)

    -- UPPER HALF
    OOUx, OOUy = rotatexy(OOx, OOy, t)
    OTUx, OTUy = rotatexy(OTx, OTy, t)

    ETUx, ETUy = rotatexy(ETx, ETy, t)
    IEUx, IEUy = rotatexy(IEx, IEy, t)

    BIUx, BIUy = rotatexy(BIx, BIy, t)
    BTUx, BTUy = rotatexy(BTx, BTy, t)

    RLUx, RLUy = rotatexy(RLx, RLy, t)
    RSUx, RSUy = rotatexy(RSx, RSy, t)

    -- Nodes
    mi_addnode(RSlx, RSly)
    mi_addnode(RLlx, RLly)
    mi_addnode(BTlx, BTly)
    mi_addnode(BIlx, BIly)
    mi_addnode(IElx, IEly)
    mi_addnode(ETlx, ETly)
    mi_addnode(OTlx, OTly)
    mi_addnode(OOlx, OOly)

    mi_addnode(OOUx, OOUy)
    mi_addnode(OTUx, OTUy)
    mi_addnode(ETUx, ETUy)
    mi_addnode(IEUx, IEUy)
    mi_addnode(BIUx, BIUy)
    mi_addnode(BTUx, BTUy)
    mi_addnode(RLUx, RLUy)
    mi_addnode(RSUx, RSUy)

    -- Lower root fillet
    mi_addarc(
        RLlx, RLly,
        RSlx, RSly,
        root_fillet_angle,
        1
    )

    -- Lower tooth body
    mi_addsegment(
        RLlx, RLly,
        BTlx, BTly
    )

    -- Lower body-to-inner-circle fillet
    mi_addarc(
        BIlx, BIly,
        BTlx, BTly,
        body_tip_fillet_angle,
        1
    )

    -- Lower concentric inner shoe arc
    mi_addarc(
        IElx, IEly,
        BIlx, BIly,
        inner_shoe_arc_angle,
        1
    )

    -- Lower inner end fillet
    mi_addarc(
        IElx, IEly,
        ETlx, ETly,
        inner_end_fillet_angle,
        1
    )

    -- Lower radial end face
    mi_addsegment(
        ETlx, ETly,
        OTlx, OTly
    )

    -- Lower outer end fillet
    mi_addarc(
        OTlx, OTly,
        OOlx, OOly,
        outer_end_fillet_angle,
        1
    )

    -- Outer concentric shoe face
    mi_addarc(
        OOlx, OOly,
        OOUx, OOUy,
        2*outer_face_half_angle,
        1
    )

    -- Upper outer end fillet
    mi_addarc(
        OOUx, OOUy,
        OTUx, OTUy,
        outer_end_fillet_angle,
        1
    )

    -- Upper radial end face
    mi_addsegment(
        OTUx, OTUy,
        ETUx, ETUy
    )

    -- Upper inner end fillet
    mi_addarc(
        ETUx, ETUy,
        IEUx, IEUy,
        inner_end_fillet_angle,
        1
    )

    -- Upper concentric inner shoe arc
    mi_addarc(
        BIUx, BIUy,
        IEUx, IEUy,
        inner_shoe_arc_angle,
        1
    )

    -- Upper body-to-inner-circle fillet
    mi_addarc(
        BTUx, BTUy,
        BIUx, BIUy,
        body_tip_fillet_angle,
        1
    )

    -- Upper tooth body
    mi_addsegment(
        BTUx, BTUy,
        RLUx, RLUy
    )

    -- Upper root fillet
    mi_addarc(
        RSUx, RSUy,
        RLUx, RLUy,
        root_fillet_angle,
        1
    )

    -- Circular slot-bottom arc to next tooth
    next_t = (i + 1) * slot_pitch
    NextRSlx, NextRSly =
        rotatexy(RSx, -RSy, next_t)

    mi_addnode(NextRSlx, NextRSly)

    mi_addarc(
        RSUx, RSUy,
        NextRSlx, NextRSly,
        slot_bottom_angle,
        1
    )
end

else
    error(
        "Unknown tooth_shoe_profile: "
        .. tooth_shoe_profile
    )
end


-- ============================================================
-- CURVED MAGNETS WITH PARALLEL SIDES + EXPLICIT CORNER FILLETS
--
-- IMPORTANT:
--   Do not use mi_createradius() here.  XFEMM's implementation can
--   trim the adjoining geometry without reliably creating the
--   replacement fillet arc.  The tangent points and fillet arcs are
--   therefore constructed explicitly.  The same geometry works in
--   both FEMM and XFEMM.
--
-- Optional XING2-style centre slot:
--   circular notch cut into the inner magnet face, with the two
--   mouth fillets also constructed explicitly as circle-circle
--   tangent fillets.
--
-- depth <= 0:
--   exact 180-degree main notch arc before mouth fillets.
--
-- depth > 0:
--   requested radial depth from the original inner face at the
--   pole centre to the deepest point of the notch.
-- ============================================================

R_mag_outer = R_mag_inner + mag_thick
half_mag = mag_width / 2

x_inner = sqrt(R_mag_inner*R_mag_inner - half_mag*half_mag)
x_outer = sqrt(R_mag_outer*R_mag_outer - half_mag*half_mag)

a_inner = asin(half_mag / R_mag_inner) / DEG
a_outer = asin(half_mag / R_mag_outer) / DEG

-- Add the minor circular arc between two points about a known centre.
-- The endpoint order is selected automatically so mi_addarc() always
-- receives the smaller positive CCW angle.
function add_minor_arc_about(
    x1, y1,
    x2, y2,
    cx, cy,
    maxseg
)
    v1x = x1 - cx
    v1y = y1 - cy
    v2x = x2 - cx
    v2y = y2 - cy

    r1 = sqrt(v1x*v1x + v1y*v1y)
    r2 = sqrt(v2x*v2x + v2y*v2y)

    if r1 <= 0 or r2 <= 0 then
        error("Cannot create arc about a zero-radius vector")
    end

    cang = (v1x*v2x + v1y*v2y) / (r1*r2)
    if cang > 1 then cang = 1 end
    if cang < -1 then cang = -1 end

    ang = acos(cang) / DEG
    cross = v1x*v2y - v1y*v2x

    if cross >= 0 then
        mi_addarc(x1, y1, x2, y2, ang, maxseg)
    else
        mi_addarc(x2, y2, x1, y1, ang, maxseg)
    end
end

-- ------------------------------------------------------------
-- Main magnet corner treatment
--
-- In local pole coordinates the magnet side faces are y=+/-half_mag.
-- The two stator/air-gap-facing corners retain the requested fillet
-- radius and are externally tangent to R_mag_inner.
--
-- The two back-iron-facing corners are intentionally SHARP.  Keeping
-- the outer corner on the exact R_mag_outer / R_back_inner junction
-- avoids the microscopic three-material PSLG slivers produced when an
-- outer fillet is tangent to the coincident back-iron inner circle.
-- ------------------------------------------------------------

if mag_corner_r < 0 then
    error("Magnet corner radius cannot be negative")
end

-- Outer corners are always sharp.
outer_corner_x = x_outer
outer_corner_y = half_mag

if mag_corner_r == 0 then
    inner_corner_circle_x = x_inner
    inner_corner_circle_y = half_mag
    inner_corner_side_x = x_inner
else
    if mag_corner_r >= half_mag then
        error("Magnet corner radius is too large for magnet width")
    end

    inner_corner_center_y = half_mag - mag_corner_r
    inner_corner_center_r = R_mag_inner + mag_corner_r

    inner_corner_center_x =
        sqrt(
            inner_corner_center_r*inner_corner_center_r
            - inner_corner_center_y*inner_corner_center_y
        )

    inner_corner_circle_x =
        inner_corner_center_x
        * R_mag_inner / inner_corner_center_r

    inner_corner_circle_y =
        inner_corner_center_y
        * R_mag_inner / inner_corner_center_r

    inner_corner_side_x = inner_corner_center_x

    if outer_corner_x <= inner_corner_side_x then
        error("Magnet inner corner fillets consume the entire side face")
    end
end

-- FEMM embeds Lua 4.0: use 1/0, not true/false.
magnet_center_slot_enabled = 0
magnet_center_slot_width = 0
magnet_center_slot_depth = 0
magnet_center_slot_fillet = 0

if magnet_center_slot_enabled == 1 then

    if magnet_center_slot_width <= 0 then
        error("Magnet centre slot width must be > 0")
    end

    if magnet_center_slot_width >= mag_width then
        error(
            "Magnet centre slot width must be smaller than magnet width"
        )
    end

    slot_half_w = magnet_center_slot_width / 2

    if slot_half_w >= R_mag_inner then
        error("Magnet centre slot is wider than its inner radius")
    end

    -- Sharp theoretical mouth points on the original inner circle.
    slot_mouth_x =
        sqrt(
            R_mag_inner*R_mag_inner
            - slot_half_w*slot_half_w
        )

    slot_mouth_angle = asin(slot_half_w / R_mag_inner) / DEG

    -- The usable inner-face arc actually ends at the explicit main
    -- magnet corner fillet tangent, not at the original sharp corner.
    inner_corner_tangent_angle =
        atan2(inner_corner_circle_y, inner_corner_circle_x) / DEG

    if slot_mouth_angle >= inner_corner_tangent_angle then
        error(
            "Magnet centre slot overlaps the magnet side/corner fillet region"
        )
    end

    if magnet_center_slot_depth <= 0 then
        -- Exact 180-degree main arc before mouth fillets.
        slot_sag = slot_half_w
        slot_deep_x = slot_mouth_x + slot_sag
        slot_effective_depth = slot_deep_x - R_mag_inner
    else
        slot_deep_x = R_mag_inner + magnet_center_slot_depth
        slot_sag = slot_deep_x - slot_mouth_x

        if slot_sag <= 0 then
            error(
                "Magnet centre slot depth does not reach inside the magnet"
            )
        end

        slot_effective_depth = magnet_center_slot_depth
    end

    if slot_deep_x >= R_mag_outer then
        error("Magnet centre slot cuts completely through the magnet")
    end

    -- Circle through lower mouth, upper mouth, and deepest point.
    slot_circle_r =
        (slot_half_w*slot_half_w + slot_sag*slot_sag)
        / (2*slot_sag)

    slot_circle_cx = slot_deep_x - slot_circle_r

    -- The unfilleted half-arc angle is retained for diagnostics.
    slot_half_arc_angle =
        2 * atan(slot_sag / slot_half_w) / DEG

    if magnet_center_slot_fillet < 0 then
        error("Magnet centre slot fillet radius cannot be negative")
    end

    if magnet_center_slot_fillet == 0 then
        slot_inner_tangent_x = slot_mouth_x
        slot_inner_tangent_y = slot_half_w
        slot_notch_tangent_x = slot_mouth_x
        slot_notch_tangent_y = slot_half_w
    else
        -- Explicit circle-circle mouth fillet.  The fillet is external
        -- to both the stator-side inner magnet circle and the notch
        -- circle, so its centre is the intersection of their radii
        -- expanded by the requested fillet radius.
        slot_fillet_ra = R_mag_inner + magnet_center_slot_fillet
        slot_fillet_rb = slot_circle_r + magnet_center_slot_fillet
        slot_centres_d = abs(slot_circle_cx)

        if slot_centres_d <= 0 then
            error("Degenerate centre-slot circle geometry")
        end

        slot_fillet_center_x =
            (
                slot_fillet_ra*slot_fillet_ra
                - slot_fillet_rb*slot_fillet_rb
                + slot_centres_d*slot_centres_d
            ) / (2*slot_centres_d)

        slot_fillet_y2 =
            slot_fillet_ra*slot_fillet_ra
            - slot_fillet_center_x*slot_fillet_center_x

        if slot_fillet_y2 <= 0 then
            error("Centre-slot mouth fillet has no valid tangent solution")
        end

        slot_fillet_center_y = sqrt(slot_fillet_y2)

        slot_inner_tangent_x =
            slot_fillet_center_x
            * R_mag_inner / slot_fillet_ra

        slot_inner_tangent_y =
            slot_fillet_center_y
            * R_mag_inner / slot_fillet_ra

        slot_notch_tangent_x =
            slot_circle_cx
            + (
                slot_fillet_center_x - slot_circle_cx
              ) * slot_circle_r / slot_fillet_rb

        slot_notch_tangent_y =
            slot_fillet_center_y
            * slot_circle_r / slot_fillet_rb

        slot_inner_tangent_angle =
            atan2(slot_inner_tangent_y, slot_inner_tangent_x) / DEG

        if slot_inner_tangent_angle >= inner_corner_tangent_angle then
            error(
                "Centre-slot mouth fillet overlaps the magnet corner fillet"
            )
        end
    end

    print("")
    print("Magnet centre slot enabled")
    print("  width                  = ", magnet_center_slot_width, " mm")
    print(
        "  requested radial depth = ",
        magnet_center_slot_depth,
        " mm (<=0 means semicircle)"
    )
    print("  effective centre depth = ", slot_effective_depth, " mm")
    print("  main arc total angle   = ", 2*slot_half_arc_angle, " deg")
    print("  mouth fillet radius    = ", magnet_center_slot_fillet, " mm")
end


for i = 0, Npoles-1 do

    t = i * pole_pitch

    -- Main-corner tangent points, local then rotated.
    ilc_x, ilc_y = rotatexy(
        inner_corner_circle_x,
        -inner_corner_circle_y,
        t
    )
    iuc_x, iuc_y = rotatexy(
        inner_corner_circle_x,
        inner_corner_circle_y,
        t
    )

    ils_x, ils_y = rotatexy(inner_corner_side_x, -half_mag, t)
    ius_x, ius_y = rotatexy(inner_corner_side_x,  half_mag, t)

    ol_x, ol_y = rotatexy(outer_corner_x, -outer_corner_y, t)
    ou_x, ou_y = rotatexy(outer_corner_x,  outer_corner_y, t)

    mi_addnode(ilc_x, ilc_y)
    mi_addnode(iuc_x, iuc_y)
    mi_addnode(ils_x, ils_y)
    mi_addnode(ius_x, ius_y)
    mi_addnode(ol_x, ol_y)
    mi_addnode(ou_x, ou_y)

    if magnet_center_slot_enabled == 1 then

        sli_x, sli_y = rotatexy(
            slot_inner_tangent_x,
            -slot_inner_tangent_y,
            t
        )
        sui_x, sui_y = rotatexy(
            slot_inner_tangent_x,
            slot_inner_tangent_y,
            t
        )

        sln_x, sln_y = rotatexy(
            slot_notch_tangent_x,
            -slot_notch_tangent_y,
            t
        )
        sun_x, sun_y = rotatexy(
            slot_notch_tangent_x,
            slot_notch_tangent_y,
            t
        )

        sd_x, sd_y = rotatexy(slot_deep_x, 0, t)
        snc_x, snc_y = rotatexy(slot_circle_cx, 0, t)

        mi_addnode(sli_x, sli_y)
        mi_addnode(sui_x, sui_y)
        mi_addnode(sln_x, sln_y)
        mi_addnode(sun_x, sun_y)
        mi_addnode(sd_x, sd_y)

        -- Original inner magnet face below the notch.
        add_minor_arc_about(
            ilc_x, ilc_y,
            sli_x, sli_y,
            0, 0,
            1
        )

        if magnet_center_slot_fillet > 0 then
            slfc_x, slfc_y = rotatexy(
                slot_fillet_center_x,
                -slot_fillet_center_y,
                t
            )

            add_minor_arc_about(
                sli_x, sli_y,
                sln_x, sln_y,
                slfc_x, slfc_y,
                1
            )
        end

        -- Main circular notch, kept as two arcs through the deepest point.
        add_minor_arc_about(
            sln_x, sln_y,
            sd_x, sd_y,
            snc_x, snc_y,
            1
        )

        add_minor_arc_about(
            sd_x, sd_y,
            sun_x, sun_y,
            snc_x, snc_y,
            1
        )

        if magnet_center_slot_fillet > 0 then
            sufc_x, sufc_y = rotatexy(
                slot_fillet_center_x,
                slot_fillet_center_y,
                t
            )

            add_minor_arc_about(
                sun_x, sun_y,
                sui_x, sui_y,
                sufc_x, sufc_y,
                1
            )
        end

        -- Original inner magnet face above the notch.
        add_minor_arc_about(
            sui_x, sui_y,
            iuc_x, iuc_y,
            0, 0,
            1
        )

    else
        add_minor_arc_about(
            ilc_x, ilc_y,
            iuc_x, iuc_y,
            0, 0,
            1
        )
    end

    -- Upper inner corner fillet.
    if mag_corner_r > 0 then
        ciu_x, ciu_y = rotatexy(
            inner_corner_center_x,
            inner_corner_center_y,
            t
        )
        add_minor_arc_about(
            iuc_x, iuc_y,
            ius_x, ius_y,
            ciu_x, ciu_y,
            1
        )
    end

    -- Upper side face ends directly at the sharp back-iron-facing corner.
    mi_addsegment(ius_x, ius_y, ou_x, ou_y)

    -- Outer curved face runs directly between the two sharp corners.
    -- It remains coincident with R_back_inner, as in the square-corner
    -- geometry that meshes reliably in FEMM/XFEMM.
    add_minor_arc_about(
        ol_x, ol_y,
        ou_x, ou_y,
        0, 0,
        1
    )

    -- Lower side face starts directly at the sharp outer corner.
    mi_addsegment(ol_x, ol_y, ils_x, ils_y)

    -- Lower inner corner fillet.
    if mag_corner_r > 0 then
        cil_x, cil_y = rotatexy(
            inner_corner_center_x,
            -inner_corner_center_y,
            t
        )
        add_minor_arc_about(
            ils_x, ils_y,
            ilc_x, ilc_y,
            cil_x, cil_y,
            1
        )
    end
end

-- ============================================================
-- ROTOR BACK IRON
-- ============================================================

addcircle(R_back_inner)
addcircle(R_back_outer)

-- ============================================================
-- OUTER AIR DOMAIN
-- ============================================================

addcircle(R_boundary)


-- ============================================================
-- SLIDING-BAND ROTOR MOTION GEOMETRY
--
-- Split the 0.15 mm physical air gap into:
--
--   stator-side air   : 1/3 gap
--   sliding band      : 1/3 gap (special unmeshed region)
--   rotor-side air    : 1/3 gap
--
-- The rotor geometry itself never has to move or remesh.
-- ============================================================

airgap = R_mag_inner - R_stator

R_band_inner =
    R_stator + airgap/3

R_band_outer =
    R_stator + 2*airgap/3

addcircle(R_band_inner)
addcircle(R_band_outer)


-- ============================================================
-- WINDING REGION GEOMETRY
--
-- Each physical slot contains one side of each of the two coils
-- on the adjacent teeth.  We divide each slot into two regions
-- with a fictitious nonmagnetic boundary.
--
-- The extra boundaries do not represent steel or insulation;
-- they merely let FEMM assign two different circuit/turn
-- properties inside one physical slot.
-- ============================================================

Nturns = 15

-- Close the mouth of each winding cavity along the inner shoe
-- tangency points, then divide the slot along its centreline.
for i = 0, Nslots-1 do

    tooth_a = i
    tooth_b = mod(i + 1, Nslots)

    ta = tooth_a * slot_pitch
    tb = tooth_b * slot_pitch

    -- Upper inner-shoe tangency on lower-angle tooth
    C1x, C1y = rotatexy(IEx,  IEy, ta)

    -- Lower inner-shoe tangency on upper-angle tooth
    C2x, C2y = rotatexy(IEx, -IEy, tb)

    -- Fictitious closure across the slot mouth
    mi_addnode(C1x, C1y)
    mi_addnode(C2x, C2y)
    mi_addsegment(C1x, C1y, C2x, C2y)

    -- Midpoint of closure
    CMx = (C1x + C2x) / 2
    CMy = (C1y + C2y) / 2
    mi_addnode(CMx, CMy)

    -- Midpoint of the circular slot-bottom arc
    slot_center_angle = (i + 0.5) * slot_pitch
    SBx, SBy = polar(R_slotbottom, slot_center_angle)

    -- Adding a node on the existing slot-bottom arc splits it.
    mi_addnode(SBx, SBy)

    -- Divider between the two coil sides
    mi_addsegment(SBx, SBy, CMx, CMy)
end


-- ============================================================
-- WINDING PHASE HELPERS
--
-- For 9 slots / 12 poles (6 pole pairs), adjacent tooth axes
-- are separated by:
--
--     6 * 40 deg = 240 electrical deg = -120 deg
--
-- An equivalent balanced phase assignment is therefore:
--
--     tooth:  0 1 2 3 4 5 6 7 8
--     phase:  A C B A C B A C B
--
-- All tooth coils use the same physical winding sense.
-- ============================================================

function tooth_phase(k)
    if k == 0 then return "A"
    elseif k == 1 then return "C"
    elseif k == 2 then return "B"
    elseif k == 3 then return "A"
    elseif k == 4 then return "C"
    elseif k == 5 then return "B"
    elseif k == 6 then return "A"
    elseif k == 7 then return "C"
    elseif k == 8 then return "B"
    end
    error("Invalid tooth index")
end

function tooth_sign(k)
    if k == 0 then return 1
    elseif k == 1 then return 1
    elseif k == 2 then return 1
    elseif k == 3 then return 1
    elseif k == 4 then return 1
    elseif k == 5 then return 1
    elseif k == 6 then return 1
    elseif k == 7 then return 1
    elseif k == 8 then return 1
    end
    error("Invalid tooth index")
end

-- ============================================================
-- CONNECTION-AWARE CURRENT HELPER
--
-- Iq_line_peak is the PEAK FUNDAMENTAL TERMINAL LINE CURRENT.
--
-- Wye:
--   winding branch currents == line currents.
--
-- Delta:
--   A/B/C circuits represent the AB / BC / CA winding branches.
--   The zero-circulating-current branch solution is:
--
--     iAB = (IA - IB)/3
--     iBC = (IB - IC)/3
--     iCA = (IC - IA)/3
--
-- which gives branch-current peak = line-current peak / sqrt(3)
-- with the correct 30-degree phase shift.
-- ============================================================

connection = "delta"

function set_line_iq_currents(Iq_line_peak, theta_mech)

    theta_e =
        (Npoles/2)
        * theta_mech

    -- theta_e is referenced to the physical phase-winding axes.
    -- For Wye, terminal line current is phase current.
    -- For Delta, branch current is shifted +30 electrical degrees
    -- relative to terminal line current.  Therefore shift the commanded
    -- terminal line-current vector -30 degrees so the resulting Delta
    -- branch-current vector lies on the same q-axis as the Wye case.
    theta_line = theta_e

    if connection == "delta" then
        theta_line = theta_e - 30
    end

    IlineA =
        Iq_line_peak
        * sin(theta_line*DEG)

    IlineB =
        Iq_line_peak
        * sin((theta_line - 120)*DEG)

    IlineC =
        Iq_line_peak
        * sin((theta_line + 120)*DEG)

    if connection == "wye" then

        IwindA = IlineA
        IwindB = IlineB
        IwindC = IlineC

    elseif connection == "delta" then

        IwindA =
            (IlineA - IlineB) / 3

        IwindB =
            (IlineB - IlineC) / 3

        IwindC =
            (IlineC - IlineA) / 3

    else
        error("Unknown winding connection: " .. connection)
    end

    mi_modifycircprop("A", 1, IwindA)
    mi_modifycircprop("B", 1, IwindB)
    mi_modifycircprop("C", 1, IwindC)

    return
        IlineA, IlineB, IlineC,
        IwindA, IwindB, IwindC
end



-- ============================================================
-- MATERIALS
-- ============================================================

-- Built-in FEMM material library entries.
-- M-19 Steel is our first-pass stator lamination approximation.
-- 1010 Steel is our first-pass rotor back-iron approximation.
mi_getmaterial("Air")
mi_getmaterial("1010 Steel")

mi_addmaterial(
    "StatorSteel",
    1, 1,
    0,
    0,
    1.9,
    0.2,
    0,
    0.93,
    0,
    0, 0,
    0, 0
)
mi_addbhpoint("StatorSteel", 0, 0)
mi_addbhpoint("StatorSteel", 0.0478945315895419, 15.120714)
mi_addbhpoint("StatorSteel", 0.0957890631790838, 22.718292)
mi_addbhpoint("StatorSteel", 0.143683594768626, 27.842733)
mi_addbhpoint("StatorSteel", 0.191578126358168, 31.871434)
mi_addbhpoint("StatorSteel", 0.239472657947709, 35.365044)
mi_addbhpoint("StatorSteel", 0.287367189537251, 38.600588)
mi_addbhpoint("StatorSteel", 0.335261721126793, 41.736202)
mi_addbhpoint("StatorSteel", 0.383156252716335, 44.873979)
mi_addbhpoint("StatorSteel", 0.431050784305877, 48.087807)
mi_addbhpoint("StatorSteel", 0.478945315895419, 51.437236)
mi_addbhpoint("StatorSteel", 0.526839847484961, 54.975221)
mi_addbhpoint("StatorSteel", 0.574734379074503, 58.752993)
mi_addbhpoint("StatorSteel", 0.622628910664044, 62.823644)
mi_addbhpoint("StatorSteel", 0.670523442253586, 67.245285)
mi_addbhpoint("StatorSteel", 0.718417973843128, 72.084407)
mi_addbhpoint("StatorSteel", 0.76631250543267, 77.4201)
mi_addbhpoint("StatorSteel", 0.814207037022212, 83.35002)
mi_addbhpoint("StatorSteel", 0.862101568611754, 89.999612)
mi_addbhpoint("StatorSteel", 0.909996100201296, 97.537353)
mi_addbhpoint("StatorSteel", 0.957890631790838, 106.201406)
mi_addbhpoint("StatorSteel", 1.00578516338038, 116.348464)
mi_addbhpoint("StatorSteel", 1.05367969496992, 128.547329)
mi_addbhpoint("StatorSteel", 1.10157422655946, 143.765431)
mi_addbhpoint("StatorSteel", 1.14946875814901, 163.754168)
mi_addbhpoint("StatorSteel", 1.19736328973855, 191.868158)
mi_addbhpoint("StatorSteel", 1.24525782132809, 234.833507)
mi_addbhpoint("StatorSteel", 1.29315235291763, 306.509769)
mi_addbhpoint("StatorSteel", 1.34104688450717, 435.255202)
mi_addbhpoint("StatorSteel", 1.38894141609671, 674.911968)
mi_addbhpoint("StatorSteel", 1.43683594768626, 1108.325569)
mi_addbhpoint("StatorSteel", 1.4847304792758, 1813.085468)
mi_addbhpoint("StatorSteel", 1.53262501086534, 2801.217421)
mi_addbhpoint("StatorSteel", 1.58051954245488, 4053.653117)
mi_addbhpoint("StatorSteel", 1.62841407404442, 5591.10689)
mi_addbhpoint("StatorSteel", 1.67630860563397, 7448.318413)
mi_addbhpoint("StatorSteel", 1.72420313722351, 9708.81567)
mi_addbhpoint("StatorSteel", 1.77209766881305, 12486.931615)
mi_addbhpoint("StatorSteel", 1.81999220040259, 16041.483644)
mi_addbhpoint("StatorSteel", 1.86788673199213, 21249.420624)
mi_addbhpoint("StatorSteel", 1.91578126358168, 31313.495878)
mi_addbhpoint("StatorSteel", 1.96367579517122, 53589.446877)
mi_addbhpoint("StatorSteel", 2.01157032676076, 88477.484601)
mi_addbhpoint("StatorSteel", 2.0594648583503, 124329.41054)
mi_addbhpoint("StatorSteel", 2.10735938993984, 159968.5693)
mi_addbhpoint("StatorSteel", 2.15525392152938, 197751.604272)
mi_addbhpoint("StatorSteel", 2.20314845311893, 234024.751347)

-- Homogenized winding material for the DC magnetostatic model.
-- mu_r = 1, like copper.  Conductivity is included for future
-- use but has no magnetic effect in this zero-frequency solve.
mi_addmaterial(
    "Winding",
    1, 1,           -- mu_x, mu_y
    0,              -- Hc
    0,              -- source J; circuit property supplies current
    58,             -- conductivity, MS/m
    0, 0,           -- lamination thickness / hysteresis
    1,              -- fill factor
    0,              -- lamination type
    0, 0,
    0, 0
)

-- Three series-connected phase circuits.
-- Initial currents are zero; the automated solves at the end
-- modify these values.
mi_addcircprop("A", 0, 1)
mi_addcircprop("B", 0, 1)
mi_addcircprop("C", 0, 1)


-- Permanent magnet generated from definition
mag_mu_r = 1.03877103145574
mag_Hc = 1087824.0360331
mi_addmaterial("PermanentMagnet", mag_mu_r, mag_mu_r, mag_Hc, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)


-- ============================================================
-- BLOCK LABELS / MATERIAL ASSIGNMENTS
-- ============================================================

-- Group numbers:
--   0 = air / fixed background
--   1 = stator
--   2 = complete rotor (back iron + magnets)
--
-- Keeping all rotor pieces in one group will be useful later
-- when we rotate the rotor programmatically.


-- ------------------------------------------------------------
-- Central bore air
-- ------------------------------------------------------------

mi_addblocklabel(0, 0)
mi_selectlabel(0, 0)

mi_setblockprop(
    "Air",
    0,          -- manual mesh size
    0.30,       -- mm
    "",
    0,
    0,
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Working air gap + spaces between magnets
--
-- The sliding-band region itself must remain unlabeled/unmeshed.
-- Put one Air label on each side of the band.
-- ------------------------------------------------------------

-- Stator-side air
air_inner_r =
    (R_stator + R_band_inner) / 2

air_inner_x, air_inner_y =
    polar(air_inner_r, 5.0)

mi_addblocklabel(air_inner_x, air_inner_y)
mi_selectlabel(air_inner_x, air_inner_y)

mi_setblockprop(
    "Air",
    0,
    0.025,
    "",
    0,
    0,
    0
)

mi_clearselected()


-- Rotor-side air
air_outer_r =
    (R_band_outer + R_mag_inner) / 2

air_outer_x, air_outer_y =
    polar(air_outer_r, 5.0)

mi_addblocklabel(air_outer_x, air_outer_y)
mi_selectlabel(air_outer_x, air_outer_y)

mi_setblockprop(
    "Air",
    0,
    0.025,
    "",
    0,
    0,
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Sliding-band interior: <No Mesh>
--
-- IMPORTANT:
-- FEMM requires the region between the two Periodic Air Gap
-- boundary circles to be explicitly marked as <No Mesh>.
-- Merely leaving this annulus unlabeled causes:
--
--   "Material properties have not been defined for all regions"
--
-- FEMM serializes this special block as a hole/no-mesh point.
-- ------------------------------------------------------------

band_label_r =
    (R_band_inner + R_band_outer) / 2

band_label_x, band_label_y =
    polar(band_label_r, 17.0)

mi_addblocklabel(
    band_label_x,
    band_label_y
)

mi_selectlabel(
    band_label_x,
    band_label_y
)

mi_setblockprop(
    "<No Mesh>",
    0,
    0,
    "",
    0,
    0,
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Winding regions
--
-- Slot i lies between tooth i and tooth i+1.
--
-- The lower-angle half of the slot is the CCW side of tooth i:
--     +Nturns
--
-- The upper-angle half is the CW side of tooth i+1:
--     -Nturns
--
-- This makes the two sides of every tooth coil carry opposite
-- axial current directions, as a real loop must.
-- ------------------------------------------------------------

winding_label_r = 4.80
winding_label_offset = 6.0

for i = 0, Nslots-1 do

    slot_center = (i + 0.5) * slot_pitch

    -- Half-slot adjacent to tooth i
    tooth_lo = i
    phase_lo = tooth_phase(tooth_lo)

    WLx, WLy =
        polar(
            winding_label_r,
            slot_center - winding_label_offset
        )

    mi_addblocklabel(WLx, WLy)
    mi_selectlabel(WLx, WLy)

    mi_setblockprop(
        "Winding",
        0,
        0.10,
        phase_lo,
        0,
        3,              -- winding group
        tooth_sign(tooth_lo) * Nturns
    )

    mi_clearselected()


    -- Half-slot adjacent to tooth i+1
    tooth_hi = mod(i + 1, Nslots)
    phase_hi = tooth_phase(tooth_hi)

    WHx, WHy =
        polar(
            winding_label_r,
            slot_center + winding_label_offset
        )

    mi_addblocklabel(WHx, WHy)
    mi_selectlabel(WHx, WHy)

    mi_setblockprop(
        "Winding",
        0,
        0.10,
        phase_hi,
        0,
        3,
        -tooth_sign(tooth_hi) * Nturns
    )

    mi_clearselected()
end


-- ------------------------------------------------------------
-- Outside air
-- ------------------------------------------------------------

outer_air_x = 15.0
outer_air_y = 0.0

mi_addblocklabel(outer_air_x, outer_air_y)
mi_selectlabel(outer_air_x, outer_air_y)

mi_setblockprop(
    "Air",
    0,
    1.0,
    "",
    0,
    0,
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Stator laminations
-- ------------------------------------------------------------

stator_label_r = (R_bore + R_slotbottom) / 2
stx, sty = polar(stator_label_r, slot_pitch/2)

mi_addblocklabel(stx, sty)
mi_selectlabel(stx, sty)

mi_setblockprop(
    "StatorSteel",
    0,
    0.10,
    "",
    0,
    1,          -- stator group
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Rotor back iron
-- ------------------------------------------------------------

back_label_r = (R_back_inner + R_back_outer) / 2
backx, backy = polar(back_label_r, pole_pitch/2)

mi_addblocklabel(backx, backy)
mi_selectlabel(backx, backy)

mi_setblockprop(
    "1010 Steel",
    0,
    0.10,
    "",
    0,
    2,          -- rotor group
    0
)

mi_clearselected()


-- ------------------------------------------------------------
-- Twelve alternating, radially magnetized 902TP magnets
--
-- Absolute polarity is arbitrary.  Magnet 0 is chosen inward,
-- magnet 1 outward, then alternating around the rotor.
-- ------------------------------------------------------------

mag_label_r = (R_mag_inner + R_mag_outer) / 2

for i = 0, Npoles-1 do

    t = i * pole_pitch
    mx, my = polar(mag_label_r, t)

    mi_addblocklabel(mx, my)
    mi_selectlabel(mx, my)

    if mod(i,2) == 0 then
        magdir = t + 180     -- radially inward
    else
        magdir = t           -- radially outward
    end

    mi_setblockprop(
        "PermanentMagnet",
        0,
        0.06,
        "",
        magdir,
        2,                   -- rotor group
        0
    )

    mi_clearselected()
end


-- ============================================================
-- OUTER MAGNETIC BOUNDARY
--
-- A=0 is placed on the 25 mm outer circle.  At ~2.7 times the
-- rotor radius this is adequate for our first qualitative model.
-- We can later compare against a larger boundary / open boundary
-- formulation if desired.
-- ============================================================

mi_addboundprop(
    "A=0",
    0, 0, 0, 0,
    0, 0, 0, 0,
    0
)

-- The outer circle was created as two 180-degree arcs.
mi_selectarcsegment(0,  R_boundary)
mi_selectarcsegment(0, -R_boundary)

mi_setarcsegmentprop(
    3,
    "A=0",
    0,
    0
)

mi_clearselected()


-- ============================================================
-- PERIODIC AIR-GAP / SLIDING-BAND BOUNDARY
--
-- FEMM uses:
--   BdryFormat 6 = Periodic Air Gap
--   BdryFormat 7 = Anti-periodic Air Gap
--
-- We model the full 360-degree machine, so use Periodic Air Gap.
-- The stator is the INNER side and the rotor is the OUTER side.
-- ============================================================

mi_addboundprop(
    "SlidingBand",
    0, 0, 0, 0,     -- A0, A1, A2, Phi
    0, 0,             -- Mu, Sigma
    0, 0,             -- c0, c1
    6                 -- Periodic Air Gap
)

-- Set these explicitly as well, so there is no ambiguity about
-- which boundary type/angles FEMM stored.
mi_modifyboundprop("SlidingBand", 9, 6)
mi_modifyboundprop("SlidingBand", 10, 0)
mi_modifyboundprop("SlidingBand", 11, 0)

-- Apply the same special boundary to both circles that bound
-- the unmeshed sliding band.
mi_selectarcsegment(0,  R_band_inner)
mi_selectarcsegment(0, -R_band_inner)
mi_selectarcsegment(0,  R_band_outer)
mi_selectarcsegment(0, -R_band_outer)

mi_setarcsegmentprop(
    0.5,
    "SlidingBand",
    0,
    0
)

mi_clearselected()


print("")
print("SlidingBand BC explicitly set to:")
print("  BdryFormat = 6 (Periodic Air Gap)")
print("  InnerAngle = 0 deg")
print("  OuterAngle = 0 deg")
print("  Band region explicitly marked <No Mesh>")

-- ============================================================
-- XFEMM WARM-START SEQUENCE BENCHMARK
-- ============================================================
points = {
    {"nl_a0", 0, 0},
    {"nl_a1", 1, 0},
}

for point_index = 1, getn(points) do
    point = points[point_index]
    point_label = point[1]
    job_angle = point[2]
    job_current = point[3]

    print("XFEMM_SWEEP_POINT_BEGIN," .. point_label .. "," .. job_angle .. "," .. job_current)
    IlineA, IlineB, IlineC,
    IwindA, IwindB, IwindC =
        set_line_iq_currents(job_current, job_angle)
    mi_modifyboundprop("SlidingBand", 11, job_angle)
    mi_saveas("job.fem")
    mi_analyze(1)
    mi_loadsolution()

    T = mo_gapintegral("SlidingBand", 0)
    IA, VA, lambdaA = mo_getcircuitproperties("A")
    IB, VB, lambdaB = mo_getcircuitproperties("B")
    IC, VC, lambdaC = mo_getcircuitproperties("C")
    nelem = mo_numelements()
    nnode = mo_numnodes()

    print("XFEMM_SWEEP_RESULT," .. point_label .. "," .. job_angle .. "," .. job_current .. "," .. T .. "," .. lambdaA .. "," .. lambdaB .. "," .. lambdaC .. "," .. nelem .. "," .. nnode)
    mo_close()
end

quit()
