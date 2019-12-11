# Under Construction

docker build -t psbap-core .

docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op generate-pdbbind-dataset


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op generate-sifts-urls

docker run -v /media/apc/DATA1/Internship/processing/sifts:/sifts -v /media/apc/DATA1/Internship/processing/datasets:/url sifts --quiet --no-clobber --retry-connrefused --waitretry=1 --read-timeout=60 --timeout=60 -t 0 -i /url/pdbbind_sifts_urls.tsv
		

docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op map-uniprot-variant-to-proteins


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op map-uniprot-variant-to-proteins


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op map-pocket-residues-to-uniprot-variants


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op generate-foldx-mutations-and-structure


docker run -it -v /media/apc/DATA1/Internship/processing/foldx:/pdb foldx RepairPDB 12

docker run -it -v /media/apc/DATA1/Internship/processing/foldx:/pdb foldx BuildModel 12


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op foldx-success-report

docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op foldx-energies-report


docker run -v /media/apc/DATA1/Internship/processing/config:/config -v /media/apc/DATA1/Internship/processing:/processing -v /media/apc/DATA1/Internship/data:/data psbap-core -op foldx-energies-html-report
