naming discussion:
    - "LFC", "Log FC", "Log 2 FC"
    - "ORF" "Orf" "ORF ID"
    - "M1 Pred log Count"
    - "M1 Predicted Count"
    - "Coord"
    - "Adj Pval", "Adj P Value"
    - "Non Insertions"


replacements:
```shell
sd '"Orf"' '"ORF"' src/**/*.py
sd 'from pytransit.universal_data import universal' 'from pytransit.interfaces import gui, cli' src/**/*.py
sd 'from pytransit.universal_data import universal' 'from pytransit.interfaces import gui, cli' src/**/*.py
sd 'universal.interface == "gui"' 'gui.is_active' src/**/*.py
sd 'universal.interface != "gui"' 'not gui.is_active' src/**/*.py
sd 'universal.interface != "cli"' 'gui.is_active' src/**/*.py
sd 'universal.interface == "cli"' 'not gui.is_active' src/**/*.py
```