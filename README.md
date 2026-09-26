# MR distortion phantom

Independent CT-referenced research analysis of MRI distortion correction in a 603-GS phantom.
This is an independently implemented workflow, not the proprietary vendor algorithm and not a clinically certified device.

## Reproduce the measurements

```bash
python -m pip install -e .
python -m mr_distortion_phantom reproduce --output results
python -m mr_distortion_phantom figures --output .
python -m unittest discover -s tests -v
```

The included 13,568 measured phantom landmarks represent 16 series, nine original protocols and seven correction pairs. Each pair uses the same 848 nodes: 91 central nodes for a six-degree-of-freedom rigid fit, 757 outer nodes for evaluation. `reproduce` independently refits the measured coordinates and checks every summary metric against the reference export (tolerance 1e-9 mm). This reproduces the landmark-analysis stage; it does not redetect images from the shared landmark table.

## Analyze locally held phantom DICOM

```bash
mr-distortion-phantom analyze --phantom-only --ct /path/ct --original /path/mr --corrected /path/corrected --ct-crop 170 333 100 360 140 375 --output results/local
```

The crop is a **dataset-specific example**, in zero-based z/y/x index bounds, not a universal setting. Select the cranial grid for each CT. Supported inputs are classic single-frame DICOM, one uniformly oriented series per directory. The fixed near-axis-aligned lattice detector is restricted to the validated phantom setup. It rejects large correspondence errors; weak or highly displaced points may be excluded by the 1.2-mm lattice-index residual filter. Compressed DICOM may require an appropriate pydicom decoder. No raw DICOM, identifiers, dates, institution names, patient images or vendor PDFs are included here.

## Interpretation

The CT is an operational reference with its own uncertainty, not physical ground truth. Nominal lattice spacing is used only for point identity; CT measurements are not forced onto an ideal grid. The primary descriptive endpoint proposed for the paper is change in P95 per protocol; this was selected after viewing the data, not preregistered. Repeated scans/repositioning are not available. Nodes are not independent scans. Nine protocols come from five scanner groups; paired inference across independent centres is not supported.

Report-derived CT/CAD distributions refer to 858/859 points and do not provide matched node IDs. Their comparison with our 757 evaluation points is descriptive. Two original protocols (P08/P09) differ by approximately 0.10 mm in mean error from the reports; parameter and registration sensitivity do not resolve this difference. No vendor-equivalence claim is made.

## Public data and examples

Protocol/scanner labels are neutral codes. Phantom coordinates are translated to separate series-centred frames; rigid-fit residuals are unchanged. Real phantom crops and purely synthetic target examples are labelled as such. These files contain no human participant data. Blinding patient data alone would not establish an ethics exemption; this repository simply does not include such data. Institutional classification remains the responsibility of the research team.

`figures/`, `tables/` and `supplement/` include derived results and captions. MIT covers original code and original derived material in this repository; third-party software retains its own license. Vendor images, reports and CAD files are not relicensed or distributed.

## Citation and status

Concept-stage research release 0.1.0. No accepted paper or DOI exists for this project. Cite the repository URL and release version until a formal archive exists. AI-assisted coding and drafting were used; scientific authorship, inspection and final approval remain with the research team. See `docs/METHODS.md` and `docs/CAPTIONS.md`.
