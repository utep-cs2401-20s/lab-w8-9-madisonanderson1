class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;


  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];
    incrCodons(inCodon);
    next = null;
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
    if(aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)){
      incrCodons(inCodon);
    }
    if(next != null){
      next.addCodon(inCodon);
    } else {
      next = new AminoAcidLL(inCodon);
    }
  }

  private void incrCodons(String codon){
    for(int i = 0; i< codons.length; i++){
      if(codons[i].equals(codon)){
        counts[i]++;
        break;
      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for(int i = 0;i < counts.length; i++){
      sum += counts[i];
    }
    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if (this == null){
      return totalDiff(inList);
    }else {
      return aminoAcidCompare(next);
    }
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    if (this == null){
      return codonDiff(inList);
    }else {
      return codonCompare(next);
    }
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    char[] aaList= new char[aminoAcidList().length];
    return getAminoAcids(aaList, this, 0);
  }

  public char[] getAminoAcids(char[] list, AminoAcidLL head, int count){
    if (head == null){
      list[count] = head.aminoAcid;
      return list;
    }else {
      list[count] = head.aminoAcid;
      return getAminoAcids(list, head.next, count++);
    }
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    int[] aaList= new int[aminoAcidList().length];
    return getAminoAcidCount(aaList, this, 0);
  }

  public int[] getAminoAcidCount(int[] list, AminoAcidLL head, int count){
    if (head == null){
      list[count] = head.aminoAcid;
      return list;
    }else {
      list[count] = head.aminoAcid;
      return getAminoAcidCount(list, head.next, count++);
    }
  }
  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    return isSorted(this);
  }

  public boolean isSorted(AminoAcidLL head){
    if(head.aminoAcid > next.aminoAcid){//not sorted
      return false;
    }
    if(head.aminoAcid > next.aminoAcid){// continue in the linked list
      return isSorted(head.next);
    }
    if (head == null){
      return true;
    }
    return false;
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    int count = 3;
    AminoAcidLL iter = new AminoAcidLL(inSequence.substring(0,3));
    while (count < inSequence.length()){
      if(AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3)) == '*'){
        break;
      }
      inSequence = inSequence.substring(count);
      iter.addCodon(inSequence.substring(0,3));
    }
    return iter;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public AminoAcidLL sort(AminoAcidLL inList){
    if(inList == null || inList.next == null){
      return inList;
    }
    AminoAcidLL middle = getMiddleNode(inList);
    AminoAcidLL secondList = middle.next;
    middle.next = null;

    return merge(sort(inList), sort(secondList));
  }

  public AminoAcidLL getMiddleNode(AminoAcidLL node){
    if(node == null){
      return null;
    }

    AminoAcidLL a = node;
    AminoAcidLL b = node.next;

    while(b != null && b.next != null){
      a = a.next;
      b = b.next;
    }
    return a;
  }

  public AminoAcidLL merge(AminoAcidLL a, AminoAcidLL b){
    AminoAcidLL temp = new AminoAcidLL();
    AminoAcidLL result = temp;

    while (a != null && b != null){
      if(a.aminoAcid < b.aminoAcid){
        temp.next = a;
        a = a.next;
      } else {
        temp.next = b;
        b = b.next;
      }
      temp = temp.next;
    }
    temp.next = (a == null) ? b:a;
    return result.next;
  }
}