# Carbon-fibre tibial prosthesis: load-line analysis

Data analysis for a school research project (TIPE, French *classes préparatoires*)
on sports tibial prostheses: how does the point where the load is applied on the
carbon blade change its stiffness and its energy return?

## Question

Compression tests were run on a carbon-fibre prosthetic blade with the load applied
at three distances from the blade edge (the "load line"), corresponding to three
blade angles. The goal is to quantify, for each setting, the force/displacement
response, the stiffness, and the share of the stored energy the blade gives back.

## Method

1. **Test bench data**: force and displacement recorded on a compression machine,
   one loading/unloading cycle per setting (`data/*/raw.txt`), plus a version with the
   zero offset corrected (`data/*/corrected.txt`).
2. **Energy balance**: numerical integration (trapezoidal rule) of the force-displacement
   curve. The loading branch gives the stored energy, the unloading branch the returned
   energy; the area between them is the hysteresis loss.
3. **Stiffness**: peak force divided by maximum displacement.
4. **Offset correction**: quadratic least-squares fit `F = a x² + b x + c` on the last
   20 raw points, solved with `(AᵀA)⁻¹AᵀY` (see `figures/least_squares_method.png`).

## Results

| Angle | Load line | Peak force | Max displacement | Stiffness | Stored | Returned | Return |
|------:|----------:|-----------:|-----------------:|----------:|-------:|---------:|-------:|
| 20° | 10 mm | 1864 N | 34.6 mm | 53.8 kN/m | 32.5 J | 29.5 J | 90.7 % |
| 17° | 8 mm | 1798 N | 42.6 mm | 42.2 kN/m | 40.6 J | 38.9 J | 95.7 % |
| 14° | 6 mm | 1701 N | 51.2 mm | 33.2 kN/m | 48.0 J | 45.7 J | 95.2 % |

A load line closer to the blade edge gives a softer blade (lower stiffness, larger
displacement) that stores more energy, and the energy return stays above 90 % in all
three cases. These are measurements on a single blade and one cycle per setting, so
they illustrate a trend rather than a statistically established result.

![Load lines compared](figures/load_lines_comparison.png)

## Run it

```bash
pip install -r requirements.txt
python src/analyze_load_lines.py          # writes figures and results/summary.csv
```

## Repository layout

```
src/analyze_load_lines.py   analysis script
data/                       test bench measurements (raw and corrected)
figures/                    figures from the project
results/                    output of the script (summary.csv, plots)
```

## Contributions

The research project was done as a team of three with Hopper (the prosthesis company), and
the physical modelling was worked out together. I wrote the original analysis code
(data processing, least-squares fit, energy and stiffness computation, plots);
this repository reorganises my three original scripts into a single, documented one.
