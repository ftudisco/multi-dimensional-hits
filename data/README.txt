This README file accompanies the dataset representing the temporal multi-layer directed citation network used in the paper

##########################################################################
## Title:	Multi-dimensional HITS: An always computable ranking 	##
##       	for temporal multi-layer directed networks		##
## Authors:	Francesca Arrigo and Francesco Tudisco			##
## Year:	2018							##
##########################################################################

Edges are stored in the file edges.txt with the format:
| nodeid | nodeid | layerid | layerid |	timeid | weight |
|   i    |    j   |    h    |    k    |    t   |    w   |

where every row corresponds to an edge with weigth w, 
from node i on layer h to node j on layer k at time t.

ID of nodes are in the file authors.txt
ID of layers are in the file journals.txt
ID of times are in the file years.txt


