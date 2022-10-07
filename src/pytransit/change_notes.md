naming discussion:
    - "LFC", "Log FC", "Log 2 FC"
    - "ORF" "Orf" "ORF ID"
    - "M1 Pred log Count"
    - "M1 Predicted Count"
    - "Coord"
    - "Adj Pval", "Adj P Value"
    - "Non Insertions"


replacements:
    `sd '"Orf"' '"ORF"' src/**/*.py`
    `sd 'from pytransit.universal_data import universal' 'from pytransit.interfaces import gui, cli' src/**/*.py`