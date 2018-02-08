# kamikaze

# Resources
[addgene blog article on CRISPR Software](http://blog.addgene.org/the-crispr-software-matchmaker-a-new-tool-for-choosing-the-best-crispr-software-for-your-needs)

[CRISPR-GA](http://54.80.152.219/)

[EENdb](http://eendb.zfgenetics.org/)
[CrisprGE](http://crdd.osdd.net/servers/crisprge/)
[WGE](http://www.sanger.ac.uk/htgt/wge/)

# Construction Logic and terminology

- Step along sequence according to supplied mutagensis (usually every codon, 3bp)
- For each codon find Nearest PAM site to the codon to be edited
- Make edit
- The minimum basepair window that includes both the PAM site and the codon edit is referred to as a "slug". There are two possible configurations for a slug depending on whether the nearest PAM is up or downstream from the codon.

|EDIT|...| PAM |     or     | PAM |...|EDIT|

# Constraints

| Component        | Size (bp)         |
|:-----------------|------------------:|
| Riboswitch       | ?                 |
| Subpool primers  |  20-25            |
| Payload          | HA-[EDIT]-[PAM]-HA|
| gRNA promoter    |       35          |
| Seed Sequence    |       20          |

Gblock synthesis limit: 230

# Ideas
Choose codon substitutions based on codon lookup tables by organism

## Mutagenesis
- Basic codon substitution (NNN -> NNN) 20 options
- Truncation
- Single Deletions (XNN, NXN, NNX)
- Double Base Deletions (XXN, XNX, NXX)
- Codon Deletion (XXX)