# kamikaze

# Resources
[addgene blog article on CRISPR Software](http://blog.addgene.org/the-crispr-software-matchmaker-a-new-tool-for-choosing-the-best-crispr-software-for-your-needs)

[CRISPR-GA](http://54.80.152.219/)

[EENdb](http://eendb.zfgenetics.org/)
[CrisprGE](http://crdd.osdd.net/servers/crisprge/)
[WGE](http://www.sanger.ac.uk/htgt/wge/)

# Construction Logic and terminology

|-SLUG-| : |TGT|...|PAM|    or    |PAM|...|TGT|
    The minimum basepair window that includes both the PAM site and the target *unedited* codon is referred to as a "slug". There are two possible configurations for a slug depending on whether the nearest PAM is up or downstream from the codon.

|-ESLUG-| : |EDT|...|PAM|    or    |PAM|...|EDT|
    The minimum basepair window that includes both the PAM site and the newly *edited* codon is referred to as a "edit slug". This part doesn't truly exist but is generated as an intermediate between Slug and Payload.

|-PAYLOAD-| : |...[HA]...|EDITED-SLUG|...[HA]...|


- Step along sequence according to supplied mutagensis (usually every codon, 3bp)
- For each codon find Nearest PAM site to the codon to be edited and generate a slug
- For each Slug make the desired edit based on mutagenesis plan
- Add Homology Arms (HA) to complete generating the payload

# Edit Cassette
# Components

| Component        | Size (bp)         |
|:-----------------|------------------:|
| Riboswitch       |  43               |
| Subpool primers  |  20-25            |
| Payload          | HA-[EDIT]-[PAM]-HA|
| gRNA promoter    |       35          |
| CRISPR RNA       |       20          |
| TRACR RNA        |       106         |

Gblock synthesis limit: 230

gRNA
- CRISPR RNA: Binds to DNA
- TRACR RNA: Binds to Enzyme
- - 

HA
- syn-pam

?Riboswitch?
?self-cleaveable-pam?

3` |     20-25      |         |     35        |     20     |                    | 5'
3` | Subpool Primer | Payload | gRNA Promoter | CRISPR RNA | ~20bp of TRACR RNA | 5'



# Ideas
Choose codon substitutions based on codon lookup tables by organism

## Mutagenesis
- Basic codon substitution (NNN -> NNN) 20 options
- Truncation
- Single Deletions (XNN, NXN, NNX)
- Double Base Deletions (XXN, XNX, NXX)
- Codon Deletion (XXX)