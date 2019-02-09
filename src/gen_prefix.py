from dna import quote, asnat

def change_repair_guide_page(n):
  # Finds marker IFPFP and then changes the following7-bit number
  #         (  <IF PF P> )<    7 >       [     0  ^  0] binary-n EOF
  return 'IIPIFFCPICPICIICIPCCCPIIC' + 'IP' + 'P' + 'P' + quote(asnat(n) + 'IIC'
