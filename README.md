# Generate interactive Brillouine zone with symmetry line from POSCAR.


### Requirement:

spglib
plotly

### usage example

For given PROCAR.gaas, run

>python ./source/getsymline.py POSCAR.gaas


Then you get
* BZ.POSCAR.gaas.html: interactive view of BZ in plotly, together with symmetry lines.
* syml.POSCAR.gaas: coordinates of symmetly line

Three qlat means three inverse primitive vectors.


###
This is on top of a little older version of original seekpath.
Thus we may need to take latest version of seekpah.
