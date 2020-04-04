import org.junit.Test;

import java.sql.SQLOutput;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
public class AminoAcidLLTester extends AminoAcidLL {

    @Test
    public void testCreateRna(){ // inital test case, i would like to see it create from RNA sequence since it is one of the must critical
        // methods in this lab because without it nothing will work
      String test = "AGG*";
      AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
      char[] expected = {'R'}; //also a bonus test of amino acid list
      assertArrayEquals(expected, var.aminoAcidList()); //passed this test
    }

    @Test
    public void testCreateRna2(){// I want to be certain that the createFromRnaSequence works because it is the most important method because if it does not work nothing else will
        String test = "AGGAAG*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        // Ran into some errors when I wanted to used the aminoacidlist method with two amino acids took me a while to fix it
        // glad I could catch this error and fixed it
        char[] expected = {'R', 'K'};
        assertArrayEquals(expected, var.aminoAcidList());
    }

    @Test
    public void testAminoAcid(){ //I wanted to ensure that since I had many problems with my last test case due to the
        //length of the string I wanted to be converted I want to test if these methods work
        String test= "AGGAAGUUC*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        char[] expected = {'R', 'K', 'F'};
        assertArrayEquals(expected, var.aminoAcidList());
    }

    @Test
    public void testAminAcidList(){ // I believe that Amino acid List is one of the most important methods so I now want to make sure it works with duplicates because I had problems with it in the past
        String test= "AAGAAGUUC*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        char[] expected = {'K', 'K', 'F'};
        assertArrayEquals(expected, var.aminoAcidList());
    }

    @Test
    public void testAminoCount(){//Inital Testing of Amino Acid Count because i did it similar to AminoAcidList so I need to see if I face the same problems so I can fix them
        String test= "AGGAAGUUC*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        var.aminoAcidCounts();
        int[] expected = {1, 1, 1};
        // just as I expected this test failed so I went back and fixed it
        assertArrayEquals(expected, var.aminoAcidCounts());
    }

    @Test
    public void testAminoCount2(){ //since I was having problems with the previous test case I want to try the same method again
        // but with a different length of inputs since that was the factor that was breaking my code
        String test= "AGGAAG*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        var.aminoAcidCounts();
        int[] expected = {1, 1};
        assertArrayEquals(expected, var.aminoAcidCounts());
    }

    @Test
    public void testAminoCount3(){// I want to test this method again although I already tested it I want to test it with this length
        // since a similar method broke with this length
        String test= "AGG*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        var.aminoAcidCounts();
        // at first this method did not work but I fixed a little mistake in my code
        int[] expected = {1};
        assertArrayEquals(expected, var.aminoAcidCounts());
    }

    @Test
    public void testAminoCount4(){//I wanted to test the method when two Amino acids were added
        String test= "AGGAAGAGG*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        var.aminoAcidCounts();
        int[] expected = {2, 1, 1};
        assertArrayEquals(expected, var.aminoAcidCounts());
    }

    @Test
    public void testSort(){ // Need to test the sorting method of the list
        String test= "AGGAAGUUC*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        var.sort(var);
        char[] expected = {'F', 'K', 'R'};
        assertArrayEquals(expected, var.aminoAcidList());
    }

    @Test
    public void testIsSorted(){
        String test= "AGGAAGUUC*";
        AminoAcidLL var = AminoAcidLL.createFromRNASequence(test);
        assertEquals(false, isSorted());
    }
}

