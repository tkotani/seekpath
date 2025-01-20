# Generate interactive Brillouine zone with symmetry lines from POSCAR.

This is applicable to any POSCAR
You can see examples at https://ecalj.sakura.ne.jp/BZgetsyml/

### Requirement:

>pip install spglib plotly

Probably, this install all the dependencies to your system.

### usage example

For given PROCAR.gaas, run

>python ./source/getsymline.py POSCAR.gaas


Then you get
* BZ.POSCAR.gaas.html: interactive view of BZ in plotly, together with symmetry lines.
* syml.POSCAR.gaas: coordinates of symmetly line



Three qlat means three inverse primitive vectors.


### conversion to POSCAR
>pip install cif2cell
>cif2cell foobar.cif -p vasp --vasp-cartesian --vasp-format=5

or you can use any other tools such as pymatgen.

### NOTE
This is on top of a little older version of original seekpath.
Thus we may need to take latest version of seekpah.
For non-timerersal case, turn off with_time_reversal=True in ./source/getsymline.py.

