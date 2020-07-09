

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author sezin
 */

public class read_SubMatch_output {

    static Map<Integer, Map<Integer, Integer>> nodeIndices = new HashMap<Integer, Map<Integer, Integer>>();  //keeps node name existing in the subMatch output and its occurence Count in specific subgraph type 
    static Map<Integer, Integer> subgraphNodecount = new HashMap<Integer, Integer>(); //keeps each subgraph type and its number of nodes 
    static Map<Integer, String> nodeIdMap = new HashMap<Integer, String>(); //keeps each nodeIndex for Submatch and its node id (String Intact etc.)
    static Map<Integer, List<String>> proteinIdsOfTestSets = new HashMap<Integer, List<String>>(); //keeps protein ids of the testSet and their indices
    static HashMap<String, String> keywordsCategories = new <String, String>  HashMap();
    static Set<String> diseaseProteins = new HashSet<String>();
    static String dataset = "";

    public static void main(String[] args) throws IOException {
        int datasetCount = 1;
        String folderName = "";
        System.out.println("Select database \n1 for NCBI,\n2 for String:\n3 for IntAct:");
        int database = (new Scanner(System.in)).nextInt();
        switch (database) {
            case 1:
                dataset = "NCBI";
                break;
            case 2:
                dataset = "STRING";
                break;
            case 3:
                dataset = "IntAct";
                break;
        }
        String method = "";
        String path="";
        folderName = path + dataset;
        final File folder = new File(folderName);
          countOccurencesInAllTypesofSubgraphs(folder, dataset);


          countDiseaseOccurencesMatrixInAllTypesofSubgraphs(folder, dataset);
      
    }

    public static void readCategoryofKeywords() throws FileNotFoundException, IOException {
        String file = "KeywordsCategory.txt";
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            String arg[] = line.split("\t");
            String keyword = arg[0];

            String category = arg[1];

            keywordsCategories.put(keyword, category);
        }
    }


   public static void countOccurencesInAllTypesofSubgraphs(final File folder, String dataset) throws IOException {
        int typeCount = 0;
        String vertexIdsMap = "VertexIndexes" + dataset + "Keywords-GRAMI.txt";
        readVertexIdsMapFile(vertexIdsMap);
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.getName().contains("subgraph")) {
                continue;
            }
            typeCount++;

            int subgraphId = Integer.parseInt(fileEntry.getName());
            System.out.println(fileEntry.getName() + "..........." + dataset);
            BufferedReader bf = new BufferedReader(new FileReader(fileEntry));
            String line = bf.readLine();
            while (line != null) {
                String args[] = line.split("\\s+");
                int parse = 0;
                while (parse < args.length) {
                    int nodeId = Integer.parseInt(args[parse]);
                    if (!nodeIndices.containsKey(nodeId)) {  //Keeps all types of Subgraphs which can have same number of nodes.As many as number of files

                        Map<Integer, Integer> temp = new HashMap<Integer, Integer>();
                        int occurenceCount = 1;
                        temp.put(subgraphId, occurenceCount);
                        nodeIndices.put(nodeId, temp);
                    } else {
                        Map<Integer, Integer> temp = nodeIndices.get(nodeId); //number of nodes in subgraph, count of occurences of the node
                        if (temp.containsKey(subgraphId)) {
                            int occurenceCount = temp.get(subgraphId) + 1;
                            temp.put(subgraphId, occurenceCount);
                        } else {
                            temp.put(subgraphId, 1);
                        }
                    }
                    parse++;
                }
                line = bf.readLine();

            }
        }
        writeNodeSubgraphCountsforAllTypes(dataset, typeCount);

    }
    
        // IT KEEPS AND CONSTRUCTS A MATRIX FOR "ALL"  DISEASE COUNTS OF SUBGRAPH 1 SUB2,....SUB163
    public static void countDiseaseOccurencesMatrixInAllTypesofSubgraphs(final File folder, String dataset) throws IOException {
        int runCount = 5;
        String vertexIdsMap = "VertexIndexes" + dataset + "Keywords-GRAMI.txt";
        readVertexIdsMapFile(vertexIdsMap);
        readIdsOfTestProteins(dataset);
        readDiseaseProteins();
        for (int j = 1; j <= runCount; j++) {
            Map<String, TreeMap<Integer, Integer>> diseaseCountofProtein = new HashMap<String, TreeMap<Integer, Integer>>(); //proteinId is the Key, set keeps the number of occurences in #allSUBGRAPHTYPES(163)

            int fileCount = 0;
            String file = "withDiseaseMatrix" + dataset + "" + j + ".txt";
        //    String file = "(forbreastCancer)withDiseaseMatrix" + dataset + "" + j + ".txt";

            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            for (final File fileEntry : folder.listFiles()) {
                Map<String, Integer> diseaseCount = new HashMap<String, Integer>(); //proteinId is the Key,  keeps the number of occurences in current subgraph

                if (fileEntry.getName().contains("subgraph")) {
                    continue;
                }
                int subgraphId = Integer.parseInt(fileEntry.getName());
                System.out.println(j + "-" + fileEntry.getName() + "-" + dataset);
                BufferedReader bf = new BufferedReader(new FileReader(fileEntry));
                String line = bf.readLine();
                while (line != null) {
                    String args[] = line.split("\\s+");
                    int parse = 0;
                    List<String> proteins = new ArrayList<String>();
                    while (parse < args.length) {
                        int submatchId = Integer.parseInt(args[parse]);
                        String proteinId = nodeIdMap.get(submatchId);
                        proteins.add(proteinId);
                        parse++;
                    }

                    parse = 0;
                    while (parse < args.length) {
                        //     
                        int submatchId = Integer.parseInt(args[parse]);
                        String proteinId = nodeIdMap.get(submatchId); //DISEASE
                        List<String> proteinIds = proteinIdsOfTestSets.get(j); //test set disease proteins are not considered !!!!!!!!!!!!!!
                        if (diseaseProteins.contains(proteinId) && !proteinIds.contains(proteinId)) { // you found a train disease protein in that subgraph 
                            Iterator it = proteins.iterator(); //ALL PROTEINS IN SUBGRAPH
                            while (it.hasNext()) {
                                String protein = (String) it.next();
                                if (!proteinId.equals(protein)) { // not equal to disease protein

                                    if (!diseaseCount.containsKey(protein)) {
                                        diseaseCount.put(protein, 1);  // except disease itself increase
                                    } else { // the protein is seen in this subgraph or in the same type but different subgraph
                                        int count = diseaseCount.get(protein);
                                        diseaseCount.put(protein, ++count);
                                    }
                                }
                            }

                        }

                        parse++;
                    }
                    line = bf.readLine();
                }

                for (Map.Entry<String, Integer> entry : diseaseCount.entrySet()) {
                    String protein = entry.getKey();
                    int count = entry.getValue();
                    if (!diseaseCountofProtein.containsKey(protein)) {
                        TreeMap<Integer, Integer> map = new <Integer, Integer>TreeMap();
                        map.put(subgraphId, count);
                        diseaseCountofProtein.put(protein, map);
                    } else {
                        TreeMap<Integer, Integer> list = diseaseCountofProtein.get(protein);
                        list.put(subgraphId, count);
                    }
                }

            }
            for (Map.Entry<String, TreeMap<Integer, Integer>> entry : diseaseCountofProtein.entrySet()) {
                TreeMap<Integer, Integer> list = entry.getValue(); //counts (163)
                bw.write(entry.getKey());
                for (Map.Entry<Integer, Integer> entryCount : list.entrySet()) {

                    bw.write(" sub:" + entryCount.getKey() + " " + entryCount.getValue());
                }
                bw.write("\n");
            }
            bw.close();
        }
    }
    
 
  
    public static void readDiseaseProteins() throws FileNotFoundException, IOException {
        /* FOR ALL DISEASES*/ String file = "Results" + dataset + "diseaseProteins.txt";

        BufferedReader bw = new BufferedReader(new FileReader(file));
        String proteinId = bw.readLine();
        while (proteinId != null) {
            diseaseProteins.add(proteinId);
            proteinId = bw.readLine();
        }
        //  System.out.println(diseaseProteins.size());
    }

    public static void readDiseaseProteinsBreast() throws FileNotFoundException, IOException {
        /* FOR BREAST CANCER*/ String file = "Results" + dataset + "breastCancerProteins.txt";

        BufferedReader bw = new BufferedReader(new FileReader(file));
        String proteinId = bw.readLine();
        while (proteinId != null) {
            diseaseProteins.add(proteinId);
            proteinId = bw.readLine();
        }
        //  System.out.println(diseaseProteins.size());
    }

    public static void readIdsOfTestProteins(String dataset) throws FileNotFoundException, IOException {
        int k = 5;
        for (int j = 1; j <= k; j++) {
            String proteinIdsOfTest = "" + dataset + "" + j + "proteinIdsofTestSet.txt";

            BufferedReader bf = new BufferedReader(new FileReader(proteinIdsOfTest));
            String line = bf.readLine();
            List<String> proteinIds = new ArrayList<String>();
            while (line != null) {
                String args[] = line.split("\\s+");
                //int index = Integer.parseInt(args[0]);
                proteinIds.add(args[1]);
                //System.out.println(index + " " + args[1]);
                line = bf.readLine();
            }
            proteinIdsOfTestSets.put(j, proteinIds);

        }
        // printMap(proteinIdsOfTestSets);
    }


    public static void writeNodeSubgraphCountsforAllTypes(String dataset, int typeCount) throws FileNotFoundException, IOException {
        File file2 = new File("NodeCountsforAllTypes" + dataset + ".txt");

        if (!file2.exists()) {
            file2.createNewFile();
        }
        FileWriter fw2 = new FileWriter(file2.getAbsoluteFile());
        BufferedWriter bw2 = new BufferedWriter(fw2);
        bw2.write("SubmatchID\tDatabaseID\t");
        for (int i = 0; i < typeCount; i++) {
            bw2.write("#" + (i + 1) + "\t");
        }
        bw2.write("\n");
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : nodeIndices.entrySet()) {
            Map<Integer, Integer> map = entry.getValue(); //#2 count, #3 count,...
            String counts = "";
            for (int i = 0; i < typeCount; i++) {
                if (map.containsKey((i + 1))) {
                    counts += "\t" + map.get((i + 1));
                } else {
                    counts += "\t" + 0;
                }
            }
            bw2.write(entry.getKey() + "\t" + nodeIdMap.get(entry.getKey()) + "" + counts + "\n");
        }
        bw2.close();
    }
 

    public static void readVertexIdsMapFile(String file) throws FileNotFoundException, IOException {
        BufferedReader bf = new BufferedReader(new FileReader(file));
        String lineSubgraph = bf.readLine();
        while (lineSubgraph != null) {
            String args[] = lineSubgraph.split("\\s+");
            String nodeId = args[2];
            int nodeIndex = Integer.parseInt(args[1]);
            nodeIdMap.put(nodeIndex, nodeId);

            lineSubgraph = bf.readLine();
        }
    }


    public static <K, V> void printMap(Map<K, V> map) {
        for (Map.Entry<K, V> entry : map.entrySet()) {
            System.out.println(entry.getKey() + "\t" + entry.getValue().toString());
            //System.out.println(entry.getKey());

        }
        System.out.println(map.size());
    }

    public static void print(Set s) {
        Iterator it = s.iterator();
        while (it.hasNext()) {
            System.out.println(it.next());
        }
        System.out.println(s.size());
    }

    public static void print(List s) {
        Iterator it = s.iterator();
        while (it.hasNext()) {
            System.out.println(it.next());
        }
        System.out.println(s.size());
    }

}
