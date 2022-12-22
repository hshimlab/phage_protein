from dna_features_viewer import (GraphicFeature, GraphicRecord, CircularGraphicRecord)
features=[GraphicFeature(start=123, end=677, strand=-1, color="#EEE8AA", label="GTP cyclohydrolase I type 1 "),
GraphicFeature(start=934, end=2433, strand=+1, color="#DCDCDC", label="xanthine dehydrogenase "),
GraphicFeature(start=2435, end=4810, strand=+1, color="#FFEFD5", label="xanthine dehydrogenase "),
GraphicFeature(start=6515, end=7321, strand=-1, color="#FF1493", label="indole-3-glycerol phosphate synthase "),
GraphicFeature(start=7331, end=8386, strand=-1, color="#4682B4", label="anthranilate phosphoribosyltransferase "),
GraphicFeature(start=11548, end=14247, strand=+1, color="#D8BFD8", label="membrane alanine aminopeptidase N "),
GraphicFeature(start=28881, end=29522, strand=+1, color="#696969", label="indoleacetamide hydrolase "),
GraphicFeature(start=30123, end=31439, strand=+1, color="#FFFAFA", label="putative terminase "),
GraphicFeature(start=36004, end=36954, strand=+1, color="#00FFFF", label="capsid protein "),
GraphicFeature(start=39907, end=41262, strand=+1, color="#696969", label="tail fibers protein "),
GraphicFeature(start=44378, end=45130, strand=+1, color="#2E8B57", label="antirepressor protein "),
GraphicFeature(start=46185, end=48617, strand=+1, color="#FF69B4", label="tail length tape-measure protein 1 "),
GraphicFeature(start=48697, end=50493, strand=+1, color="#DC143C", label="tail length tape-measure protein 1 "),
GraphicFeature(start=55702, end=56289, strand=+1, color="#FFEBCD", label="secretion activator protein "),
GraphicFeature(start=58361, end=58645, strand=-1, color="#808080", label="lesion bypass DNA polymerase V "),
GraphicFeature(start=58690, end=59886, strand=-1, color="#808080", label="integrase "),
GraphicFeature(start=60374, end=62176, strand=-1, color="#00FFFF", label="Xaa-Pro aminopeptidase "),
]
record = GraphicRecord(sequence_length=65000, features=features)
ax, _ = record.plot(figure_width=30)
ax.figure.savefig("graphic_record_defined_by_hand_without HP.png")
record.plot(figure_width=30)

circular_rec = CircularGraphicRecord(sequence_length=65000, features=features)
ax2, _ = circular_rec.plot(figure_width=10)
ax2.figure.tight_layout()
ax2.figure.savefig(
    "graphic_record_defined_by_hand_circular_without HP.png", bbox_inches="tight"
)