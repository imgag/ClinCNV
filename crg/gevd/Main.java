package crg.gevd;

import javafx.util.Pair;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.io.*;
import java.util.Map;
import java.util.HashMap;


class Solver {
    private Double probOfAdditional = 0.01;
    private Double probOfMain = 1 - probOfAdditional;
    private Double probOfHomDelSNP = 0.001;
    private Double probOfMistakenglyDetectedSNP = 0.999;

    private String fileNameIn;
    private ArrayList<String> fileNamesOut;
    ArrayList<ArrayList<Double>> probs;
    ArrayList<ArrayList<Double>> negProbs;

    ArrayList<ArrayList<Double>> additionalProbs;
    ArrayList<ArrayList<Double>> weightsOfClusters;
    Integer maxNumberOfStates = 7;
    String pathToNamesConvention = "/users/so/gdemidov/CNV/snps/merged_mappings.txt";


    Double[] cn0 = {Math.log(1.0)};
    Double[] cn1 = {Math.log(1.0)};
    Double[] cn2 = {Math.log(0.4), Math.log(0.6)};
    Double[] cn3 = {Math.log(0.325), Math.log(0.295), Math.log(0.38)};
    Double[] cn4 = {Math.log(0.345), Math.log(0.22), Math.log(0.095), Math.log(0.34)};
    Double[] cn5 = {Math.log(0.39), Math.log(0.29), Math.log(0.045), Math.log(0.045), Math.log(0.23)};
    Double[] cn6 = {Math.log(0.4), Math.log(0.2), Math.log(0.03), Math.log(0.18), Math.log(0.03), Math.log(0.16)};
    ArrayList<Double[]> weights = new ArrayList<Double[]>();

    void setWeightsOfClusters() {
        weights.add(cn0);
        weights.add(cn1);
        weights.add(cn2);
        weights.add(cn3);
        weights.add(cn4);
        weights.add(cn5);
        weights.add(cn6);
    }

    static Double lbinomial(final int N, final int K) {
        Double ret = 0.0;
        for (int k = 0; k < K; k++) {
            ret += Math.log(N - k);
            ret -= Math.log(k + 1);
        }
        return ret;
    }


    private Double calculateLikelihoodBinom(ArrayList<Pair<Integer, Integer>> snpsSamplePos, Integer i) {
        Double loglik = 0.0;

        for (Pair<Integer, Integer> elem : snpsSamplePos) {
            Double likelihoodOnePosition = 0.0;
            //Integer x = elem.getKey();
            //Integer y = elem.getValue();

            for (int j = 0; j < probs.get(i).size(); j++) {
                //Double prob = probs.get(i).get(j);
                likelihoodOnePosition += Math.exp(weights.get(i)[j]  + elem.getKey() * probs.get(i).get(j) + elem.getValue() * (negProbs.get(i).get(j)));
                //likelihoodOnePosition += weights.get(i)[j] * (Math.pow(probs.get(i).get(j), elem.getKey())) * (Math.pow(1 - probs.get(i).get(j), elem.getValue()));
            }
            /*for (Double prob : additionalProbs) {
                likelihoodOnePosition += probOfAddCurrent * Math.exp(
                        lbinomial(minValue + maxValue, minValue) +
                                minValue * Math.log(prob) + maxValue * Math.log(1 - prob)
                );
            }*/
            loglik += Math.log(likelihoodOnePosition);
        }

        return 20 * loglik;
    }

    private void generateProbs() {
        probs = new ArrayList<ArrayList<Double>>();
        negProbs = new ArrayList<ArrayList<Double>>();
        ArrayList<Integer> states = new ArrayList<Integer>();
        for (int i = 0; i < maxNumberOfStates; i++){
            states.add(i);
        }
        //ArrayList<Double> probsDel = new ArrayList<Double>(1);
        //probsDel.add(Math.log(probOfHomDelSNP));
        //probs.add(0,probsDel);
        //probs.add(1,probsDel);
        for (int j = 0; j < maxNumberOfStates; j++) {
            ArrayList<Double> currentProbs = new ArrayList<Double>();
            ArrayList<Double> currentNegProbs = new ArrayList<Double>();

            //currentProbs.add(probOfMistakenglyDetectedSNP);
            currentProbs.add(Math.log(probOfHomDelSNP));
            currentNegProbs.add((Math.log(1 - probOfHomDelSNP)));
            for (int i = 1; i < states.get(j); i++) {
                Double snpAlleleBalance = (1.0 * i) / states.get(j);
                if (j > 2) {
                    currentProbs.add(Math.log(snpAlleleBalance));
                    currentNegProbs.add(Math.log(1 - snpAlleleBalance));
                }
                if (j == 2) {
                    currentProbs.add(Math.log(snpAlleleBalance));
                    currentNegProbs.add(Math.log(1 - snpAlleleBalance));
                }
            }
            probs.add(currentProbs);
            negProbs.add(currentNegProbs);
        }
        for (int i = 0; i < probs.size(); i++) {
            System.out.println(probs.get(i));
        }
    }

    private void generateAdditionalProbs() {
        additionalProbs = new ArrayList<ArrayList<Double>>();
        additionalProbs.add(new ArrayList<Double>());
        additionalProbs.get(0).add(probOfHomDelSNP);
        for (int i = 1; i < maxNumberOfStates; i++) {
            ArrayList<Double> addAlleleFreqs = new ArrayList<Double>();
            addAlleleFreqs.addAll(probs.get(i - 1));
            if (i < maxNumberOfStates - 1) {
                addAlleleFreqs.addAll(probs.get(i + 1));
            }
            if (i > 2) {
                addAlleleFreqs.add(new Double(0.01));
            }
            additionalProbs.add(addAlleleFreqs);
        }
        /*for (int i = 0; i < additionalProbs.size(); i++) {
            System.out.println(probs.get(i));
            System.out.println(additionalProbs.get(i));
            System.out.println();
        }*/

    }


    private ArrayList<Pair<Integer, Integer>> parseSnpsFromSemicolonString(String forSampleSNPs) {
        ArrayList<Pair<Integer, Integer>> snpsForSample = new ArrayList<Pair<Integer, Integer>>();
        String[] separateAlleleCounter = forSampleSNPs.split(";");
        for (String str : separateAlleleCounter) {
            String[] separateSNP = str.split(",");
            Pair<Integer, Integer> value = new Pair(Integer.parseInt(separateSNP[0]), Integer.parseInt(separateSNP[1]));
            if ((value.getKey() + value.getValue()) > 19 & (value.getKey() + value.getValue() < 150)) {
                snpsForSample.add(value);
            }
        }
        return snpsForSample;

    }
    private ArrayList<String> transformSamplesToOtherEncoding(ArrayList<String> samples) {
        ArrayList<String> transformedSamples = new ArrayList<>();
        Map<String, String> dictSampleNames = new HashMap<String, String>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(pathToNamesConvention));
            for (String line; (line = br.readLine()) != null; ) {
                String[] splittedLine = line.split("\\t");
                dictSampleNames.put(splittedLine[0], splittedLine[1]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        String topStringSampleNames = "pos\t";
        for (int i = 0; i < samples.size(); i++) {
            topStringSampleNames += dictSampleNames.get(samples.get(i));
            topStringSampleNames +="";
        }
        return transformedSamples;
    }

    public Solver(String fni, ArrayList<String> fno) {
        fileNameIn = fni;
        fileNamesOut = fno;
        Integer x = 10;
        Integer y = 5;
        Double pr = 0.33;
        Double val = lbinomial(x + y, Math.min(x,y)) + Math.min(x,y) * Math.log(pr) + Math.max(x,y) * Math.log(1 - pr);
        System.out.println("Binom likelihood is equal to...");
        System.out.println(val);
        //System.out.println(binomial(20, 5));
        //System.out.println(Math.round(Math.exp(lbinomial(20,5))));
        parseFiles();
    }

    private void parseFiles() {
        generateProbs();
        generateAdditionalProbs();
        setWeightsOfClusters();
        Integer counter = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileNameIn));
            String lineSamples = br.readLine();
            String[] sample = lineSamples.split("\t");
            ArrayList<PrintWriter> writers = new ArrayList<PrintWriter>();
            for (int i = 0; i < fileNamesOut.size(); i++) {
                FileWriter fw = new FileWriter(fileNamesOut.get(i), true);
                BufferedWriter bw = new BufferedWriter(fw);
                PrintWriter out = new PrintWriter(bw);
                writers.add(out);
            }
            for (int j = 0; j < maxNumberOfStates; j++) {
                writers.get(j).println(lineSamples);
            }
            for (String line; (line = br.readLine()) != null; ) {
                counter += 1;
                if (counter % 10000 == 0)
                    System.out.println(counter);
                String[] splittedLine = line.split("\\t");
                ArrayList<String> answerLines = new ArrayList<String>();
                for (int j = 0; j < maxNumberOfStates; j++){
                    answerLines.add(splittedLine[0] + "\t");
                }

                for (int i = 1; i < splittedLine.length; i++) {
                    if (!splittedLine[i].equals("0")) {
                        ArrayList<Pair<Integer, Integer>> snpsSamplePos = parseSnpsFromSemicolonString(splittedLine[i]);
                        if (snpsSamplePos.size() > 0 ) {
                            for (int j = 0; j < maxNumberOfStates; j++) {
                                //String shortNumber = String.format("%.2f", calculateLikelihoodBinom(snpsSamplePos, probs.get(j), additionalProbs.get(j)));
                                Integer shortNumber = calculateLikelihoodBinom(snpsSamplePos, j).intValue();
                                answerLines.set(j, answerLines.get(j) + shortNumber + "\t");
                            }
                        } else {
                            for (int j = 0; j < maxNumberOfStates; j++) {
                                answerLines.set(j, answerLines.get(j) + 0 + "\t");
                            }
                        }

                    } else {
                        for (int j = 0; j < maxNumberOfStates; j++) {
                            answerLines.set(j, answerLines.get(j) + 0 + "\t");
                        }
                    }
                }
                for (int j = 0; j < maxNumberOfStates; j++) {
                    writers.get(j).println(answerLines.get(j).substring(0, answerLines.get(j).length() - 1));
                }
            }
        } catch (FileNotFoundException e) {
            System.err.println("\nSorry, but file was not found.");
            System.err.println(e.getMessage());
            System.exit(1);
        } catch (IOException e) {
            System.err.println("\nSorry, but we had some problems with reading from your file.");
            System.err.println(e.getMessage());
            System.err.println("Please, be sure that this tool has access to your input files.");
            System.exit(1);
        }
    }
}


public class Main {

    public static void main(String[] args) {
        // write your code here
        Integer initChrom = 8;
        Integer lastChrom = 22;
        for (int i = initChrom; i <= lastChrom; i++) {
            System.out.println("We started with chromosome..." + i + "\n");
            Integer chrom = i;
            String path_in = "/users/so/gdemidov/CNV/snps/";
            String path_out = "/users/so/gdemidov/CNV/snps/likeliks/";
            String fileNameIn = path_in + "snps_het" + chrom + ".txt";
            String fileNameOut = path_out + "snps_likeliks_" + chrom;
            ArrayList<String> fileNamesOut = new ArrayList<String>(7);
            for (int j = 0; j < 7; j++) {
                fileNamesOut.add(fileNameOut + "_" + (j) + ".txt");
            }
            System.out.println(fileNameIn + "\n" + fileNamesOut);
            Solver execute = new Solver(fileNameIn, fileNamesOut);
        }
    }
}
