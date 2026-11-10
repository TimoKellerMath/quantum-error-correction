# Numbering of the group elements in the board files of four parents: the board numbers the elements of
# SmallGroup(N, id) in the order of GAP's Elements(G). For each group this prints the orders of
# IndependentGeneratorsOfAbelianGroup(G) and, in the order of Elements(G), the exponent vector of every element
# with respect to these generators, one group per line.   gap -q src/gap_elements.g | tr -d '\r' > logs/gap_elements.log
elementsLine := function(N, id)
  local G, gens;
  G := SmallGroup(N, id);
  gens := IndependentGeneratorsOfAbelianGroup(G);
  return Concatenation("SmallGroup(", String(N), ",", String(id), ") ", StructureDescription(G), " orders ",
                       ReplacedString(String(List(gens, Order)), " ", ""), " elements ",
                       ReplacedString(String(List(Elements(G), g -> IndependentGeneratorExponents(G, g))), " ", ""));
end;
SetPrintFormattingStatus("*stdout*", false);
for p in [[60, 13], [32, 21], [48, 20], [33, 1]] do
  Print(elementsLine(p[1], p[2]), "\n");
od;
QUIT;
