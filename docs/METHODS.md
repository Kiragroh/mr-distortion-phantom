# Methods and data dictionary

## Scope
Physical 603-GS MR distortion phantom, one measured planning CT, nine original MRI protocols and seven eligible corrected counterparts. Five scanner groups are represented; these are not five independent centres. No patient images, target contours, clinical plans or human biological material enter this release. Protocol labels describe source acquisition variants; scanner correction settings and exact software versions require confirmation before manuscript submission. P05 and P06 have no eligible matched corrected image in the analysed collection and are not used in paired endpoints.

## Independent raw-image implementation
Classic single-frame DICOM slices are sorted using the orientation normal. The slice vector is obtained from first-to-last image position divided by the number of intervals, preserving any constant in-plane component. Uniform orientation and slice-position residuals <0.002 mm are required. Pixel spacing, rescale slope/intercept and row/column direction vectors define physical LPS coordinates. No report coordinates enter this stage.

The cranial CT grid is cropped with explicit zero-based, half-open z/y/x bounds. In this dataset the crop was [170:333,100:360,140:375]. CT values are clipped to [−50,400]; MRI intensity is negated. Difference-of-Gaussians scales are 1.4 and 3.5 mm. Candidate maxima use an approximately 7-mm neighbourhood, boundary erosion, a response threshold at the 98th percentile, and an MRI background mask. A separable parabolic peak fit yields subvoxel locations.

Near-axis-aligned grid neighbours estimate a lattice basis. Integer indices and phase are iteratively fitted; candidates with index-fit residual ≥1.2 mm are excluded. Duplicate indices are excluded. This selection can remove strongly distorted or weak landmarks and limits extrapolation to more distorted data. Index shifts are chosen by largest common set, then rigid-fit RMS; the implementation is specific to this phantom geometry, not a universal detector. Nominal spacings (11,10.5,10 mm) define node identities and nominal radial positions. Measured CT coordinates are never replaced by an ideal lattice. The manufacturer CAD file is unavailable.

## Registration and endpoint
The centre is the componentwise median of all unique CT grid indices. Common CT/original/corrected nodes are selected within each protocol pair. In the present collection there are 848 nodes per series. A nominal radius ≤30 mm selects 91 central registration nodes. Original and corrected MRIs each receive their own six-degree-of-freedom, proper rigid least-squares fit to these CT nodes. Scale, shear and non-rigid registration are not allowed. The remaining 757 nodes form the separate evaluation set. Residuals therefore quantify deformation relative to a central rigid alignment; they are not absolute treatment-coordinate errors. Central residual means are listed separately in Table S1.

The proposed primary descriptive endpoint, selected post hoc, is ΔP95 = P95(corrected residual norm) − P95(original residual norm). It is not P95 of paired nodewise changes. NumPy percentiles use the default linear interpolation (method='linear'). Median, quartiles, P5, maximum and spatial patterns are secondary descriptions. Seven paired protocol effects are displayed separately. Spatial nodes are not independent acquisitions; no pooled node-level hypothesis test or acquisition-level confidence interval is calculated. No repeatability estimates can be obtained from this collection.

## Report comparison
Thirty-two report distributions were extracted from vector plots: 16 MRI series, each with CT and CAD references. CT reports contain 858 points and CAD reports 859. Counts, radius-shell summaries, maxima and top-100 summaries were checked against printed reports; the extraction/rounding checks differed by at most 0.00051 mm. Node identities and exact internal registration conventions are unavailable. Only aggregate descriptive comparison is justified. Similar CT/CAD distributions do not establish reference equivalence or interchangeability. Public report data contain derived numerical measurements, not the PDF reports or their layout.

## Technical verification
Four executable tests cover rigid-transform recovery, preservation of a known held-out local displacement, analytic equal-sphere overlap boundaries, and refitting all 16 series from the released coordinates. The complete raw DICOM CLI was additionally run on P02 and compared with the internal result at 1e-9-mm tolerance. This is a port/regression check, not an independent accuracy validation.

Nine original MRI volumes underwent known (+0.35,−0.45,+0.25)-mm translations with cubic resampling and redetection. These archived experiments yield P95 errors 0.081–0.093 mm but individual maxima up to 1.789 mm. Their different detected-node counts (861–863) are not the main 757-node evaluation set. These tests characterize numerical equivariance including correspondence behaviour; they do not establish absolute localization accuracy. Raw images are not bundled, so these archived image-test results cannot be rerun from the landmark CSV alone. P08 filter sensitivity at σ=1.2 and 1.6 mm retains a reduction in mean residual, but does not resolve the report discrepancy. A known non-rigid voxel-deformation experiment and repeat/repositioned scans are proposed next, not completed evidence.

## Public files and limitations
- `landmarks.csv`: 13,568 rows; protocol, state, integer grid identity i/j/k, registration/evaluation role, nominal radius, and measured CT/MR LPS coordinates in mm. Each series' coordinate set is centred independently by translation. Rigid-fit residuals are invariant to this centring.
- `protocols.json`: neutral protocol/scanner labels, manufacturer, field, DICOM TR/TE, bandwidth and voxel dimensions. DICOM TR semantics may differ across sequence families.
- `report_distributions.json`: report-derived radial residual magnitudes and summaries without locations, dates, UIDs or matched node identities.
- `technical_validation.json`: archived translation and detector-sensitivity results, not executable raw image fixtures.
- `expected_metrics.json`: frozen numeric regression reference. It is used only to check the refitted outputs.
- Static phantom images/GIFs are conventional resamplings with fixed CT markers; synthetic coverage is clearly separate and has no dose interpretation.

The CT is an operational reference, not independently calibrated physical ground truth. This single-phantom, single-reference dataset has no repeat acquisitions or uncertainty budget. Scanner-/patient-/clinical-outcome generalization and vendor equivalence are not demonstrated.
