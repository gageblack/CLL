#!/usr/bin/env python3

import sys
import os
import openpyxl
import pandas as pd

### Looks for mutated genes that show up multiple times in different group ###

## The following function takes a parsed VCF file and collects the info for each variant, 
# and adds the gene and variant to a dictionary


def get_info(patient_list, high_impact = False, onco_only = False, cll_only = True, only_AF_change = False, no_cll = False):
    if cll_only:
        onco_only = False
        no_cll = False
    gene_dict = dict()
    var_dict = dict()
    patient_genes = dict()
    genelist = set()
    for filename in os.listdir("/Users/GBlack/Desktop/CLL/2.Data/Parsed_VCFs/somatics_all/allsep"):
        patient = filename.split(".")[0]
        if filename.endswith("parsedVCF.all.tsv") and patient in patient_list:
            infile = open("/Users/GBlack/Desktop/CLL/2.Data/Parsed_VCFs/somatics_all/allsep/"+filename, "r")
            header = next(infile).strip()
            for line in infile:
                fields = line.strip().split("\t")
                gene = fields[2]
                variant = fields[3]
                var_type = fields[4]
                impact = fields[5]
                germline_DP = fields[6]
                germline_AF = fields[9]
                baseline_DP = fields[10]
                baseline_AF = fields[13]
                final_DP = fields[-4]
                final_AF = fields[-1]
                
                ### Calculate the largest change in AF for each variant
                highest = float(baseline_AF)
                lowest = float(baseline_AF)
                curr_itr = 17
                while True:
                    if float(fields[curr_itr]) > highest:
                        highest = float(fields[curr_itr])
                    if float(fields[curr_itr]) < lowest:
                        lowest = float(fields[curr_itr])
                    curr_itr = curr_itr+4
                    if curr_itr > len(fields):
                        break
                #del_AF = round(highest - lowest, 3) #Delta (change) in AF
                del_AF = round(float(final_AF)-float(baseline_AF),3)

                ## The following are filters that can be added to specify the genes you want.
                ## When calling the function you can specify which filters you want.

                # Filter low quality variants based on germline allele frequency and depth at three timepoints.
                if float(germline_AF) > 0.0095 or int(germline_DP) < 75 or int(baseline_DP) < 75 or int(final_DP) < 75:
                    continue
                
                # Only include variants that had a change of allele frequency >= 0.1 during treatment
                if only_AF_change:
                    if abs(del_AF)<0.1:
                        continue

                # Only keep variants that have a High or Moderate impact annotation.
                if high_impact and (impact == "MODIFIER" or impact == "LOW"):
                    continue

                # Only include genes that have been implicated as oncogene or suppressor.
                if onco_only:
                    onco_file = open("/Users/GBlack/Desktop/CLL/2.Data/allCancerGenes_list.txt", "r") #Can change all to just known or candidate. all include both
                    oncogenes = onco_file.read().split("\n")
                    onco_file.close()
                    if gene not in oncogenes:
                        continue
                # Only include CLL relevant genes.
                if cll_only or no_cll:
                    onco_file = open("/Users/GBlack/Desktop/CLL/2.Data/CLL_genes_updated.txt", "r")
                    oncogenes = onco_file.read().split("\n")
                    onco_file.close()
                    if cll_only and gene not in oncogenes:
                        continue
                    if no_cll and gene in oncogenes:
                        continue

                genelist.add(gene)
                
                # Add the variant to the variant and gene dictionaries, 
                # using the var or gene as the key, patient with the 
                # change in AF as the value.
                if variant in var_dict.keys():
                    var_dict[variant].add(patient)
                else: var_dict[variant] = {gene, patient}

                if gene in gene_dict.keys():
                    gene_dict[gene].add(patient+":"+str(del_AF))
                else: gene_dict[gene] = {patient+":"+str(del_AF)}

                if patient in patient_genes.keys():
                    patient_genes[patient].add(gene+":"+variant+":"+str(del_AF))
                else: patient_genes[patient] = {gene+":"+variant+":"+str(del_AF)}

    if choice == "evolution" or choice == "no_evolution":
        return (gene_dict, var_dict, patient_genes)
    if choice == "genes":
        return (gene_dict, var_dict, genelist)
    return (gene_dict, var_dict)

def print_result(dict_1, dict_2, num_needed):
    for key in sorted(dict_1.keys()):
        if key not in dict_2 and len(dict_1[key]) >= num_needed:
            print (key+": " + str(sorted(dict_1[key])))

## Read in clinical data as a dataframe
df = pd.read_excel('/Users/GBlack/Desktop/CLL/CLL_Patient_Data.xlsx', engine='openpyxl')
# Get a patient's evolution status.
#print(df.loc[df['patient'] == 252336]['evolution_status'].iloc[0])
choice = "enter"
while choice == "enter":

    print("\nWhich comparison would you like to do?")
    print("1. All Patients\n2. Responders vs non responders\n3. Evolution vs no evolution\n\
4. Response within ibrutinib patients\n5. Evolution within ibrutinib patients\n\
6. Response within acalabrutinib patients\n7. Evolution within acalabrutinib patients\n\
split: Splits into each evolution type and response\ngroups: Output in Response/evolution groups\n\
evolution: Display each variants by patient and writes it to a file \ngenes: Get all mutated genes \
    and write them to a file\nrelapse: Get relapse group")
    choice = input("Enter option: ")

    # Find somatically mutated genes found in all samples
    if choice == "1": 
        print("\nComparing all patients")
        patient_list = list(map(str,list(df['patient'])))
        gene_result, var_result = get_info(patient_list)

        print("\nGenes found in multiple patients:")
        print_result(gene_result, dict(), 1)
        print("\nSpecific variants found in multiple patients:")
        print_result(var_result, dict(), 3) #One more because the gene name is included in the list

    elif choice == 'specific':
        gene = input("Which gene do you want to get? ")
        patient_list = list(map(str,list(df['patient'])))
        gene_result, var_result = get_info(patient_list)
        #print(var_result)
        for var in var_result:
            if gene in var_result[var]:
                print(var+": "+str(var_result[var]))


    elif choice == "genes": 
        print("\ngetting all genes")
        patient_list = list(map(str,list(df['patient'])))
        gene_result, var_result, genelist = get_info(patient_list)
        outfile = open("genes.txt", "w")
        for gene in genelist:
            outfile.write(gene+'\n') 
        outfile.close() 

    elif choice == "evolution": # Displays variants in each patient
        print("Getting the genes of each patient")
        patient_list = list(map(str,list(df['patient'])))
        gene_result, var_result, pat_dict = get_info(patient_list)
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))
        #responder_list = list(map(str,list(df.loc[df['response'] == "responder"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        #relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "evolution"]['patient'])))
        outfile = open("patient.genes.txt", "w")
        print('\033[1m'+'\033[4m'+"Responders"+'\033[0m')
        for key in pat_dict.keys():
            #outfile.write(key+": " + str(sorted(pat_dict[key]))+'\n')
            if key in responder_list:# and key in evolution_list:
                print(key+": " + str(sorted(pat_dict[key])))
        print('\033[1m'+'\033[4m'+"Relapsers"+'\033[0m')
        for key in pat_dict.keys():
            #outfile.write(key+": " + str(sorted(pat_dict[key]))+'\n')
            if key in relapse_list:# and key in evolution_list:
                print(key+": " + str(sorted(pat_dict[key])))

        outfile.close()

    elif choice == "no_evolution": # Displays variants in each patient
        print("Getting the genes of each patient")
        patient_list = list(map(str,list(df['patient'])))
        gene_result, var_result, pat_dict = get_info(patient_list)
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))
        #responder_list = list(map(str,list(df.loc[df['response'] == "responder"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        #relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"]['patient'])))
        print_result(gene_result, dict(), 2)
        print('\033[1m'+'\033[4m'+"Responders"+'\033[0m')
        for key in pat_dict.keys():
            if key in responder_list and key in evolution_list:
                print(key+": " + str(sorted(pat_dict[key])))
        print('\033[1m'+'\033[4m'+"Relapsers"+'\033[0m')
        for key in pat_dict.keys():
            if key in relapse_list and key in evolution_list:
                print(key+": " + str(sorted(pat_dict[key])))
    
    elif choice == "relapse": # just get info for relapse patients
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))
        gene_result, var_result = get_info(relapse_list)
        print("\nGenes found in multiple relapse:")
        print_result(gene_result, dict(), 1)
        print_result(var_result, dict(), 3)
        
    elif choice == "responder": # just get info for relapse patients
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        gene_result, var_result = get_info(responder_list)
        print("\nGenes found in multiple responder:")
        print_result(gene_result, dict(), 1)
        print_result(var_result, dict(), 3)

    elif choice == "groups": # split output into response/evolution groups
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "evolution"]['patient'])))
        no_evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"]['patient'])))
        
        evolution_responder = list(set(responder_list) & set(evolution_list))
        evolution_relapse = list(set(relapse_list) & set(evolution_list))
        no_evolution_responder = list(set(responder_list) & set(no_evolution_list))
        no_evolution_relapse = list(set(relapse_list) & set(no_evolution_list))

        gene_result, var_result = get_info(evolution_relapse)
        print("\nGenes found in multiple evolution/relapse:")
        print_result(gene_result, dict(), 1)
        #print("\nSpecific variants found in multiple evolution/relapse:")
        #print_result(var_result, dict(), 3)

        gene_result, var_result = get_info(evolution_responder)
        print("\nGenes found in multiple evolution/responder:")
        print_result(gene_result, dict(), 1)
        #print("\nSpecific variants found in multiple evolution/responder:")
        #print_result(var_result, dict(), 3)
        
        gene_result, var_result = get_info(no_evolution_relapse)
        print("\nGenes found in multiple no_evolution/relapse:")
        print_result(gene_result, dict(), 1)
        #print("\nSpecific variants found in multiple no_evolution/relapse:")
        print_result(var_result, dict(), 3)

        gene_result, var_result = get_info(no_evolution_responder)
        print("\nGenes found in multiple no_evolution/responder:")
        print_result(gene_result, dict(), 1)
        #print("\nSpecific variants found in multiple evolution/responder:")
        print_result(var_result, dict(), 3)

    elif choice == "split": # Splits into each evolution type and response
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))
        no_evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"]['patient'])))
        replacement_list = list(map(str,list(df.loc[df['evolution_type'] == "replacement"]['patient'])))
        selection_list = list(map(str,list(df.loc[df['evolution_type'] == "selection"]['patient'])))
        evolution_list = list(map(str,list(df.loc[df['evolution_type'] == "new_clone_emergence"]['patient'])))#

        evolution_responder = list(set(responder_list) & set(evolution_list))
        evolution_relapse = list(set(relapse_list) & set(evolution_list))
        no_evolution_responder = list(set(responder_list) & set(no_evolution_list))
        no_evolution_relapse = list(set(relapse_list) & set(no_evolution_list))
        selection_responder = list(set(responder_list) & set(selection_list))
        selection_relapse = list(set(relapse_list) & set(selection_list))
        replacement_responder = list(set(responder_list) & set(replacement_list))
        replacement_relapse = list(set(relapse_list) & set(replacement_list))

        gene_result, var_result = get_info(evolution_relapse)
        print("\nGenes found in multiple evolution/relapse:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(evolution_responder)
        print("\nGenes found in multiple evolution/responder:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(no_evolution_relapse)
        print("\nGenes found in multiple no_evolution/relapse:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(no_evolution_responder)
        print("\nGenes found in multiple no_evolution/responder:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(selection_relapse)
        print("\nGenes found in multiple selection/relapse:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(selection_responder)
        print("\nGenes found in multiple selection/responder:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(replacement_relapse)
        print("\nGenes found in multiple replacement/relapse:")
        print_result(gene_result, dict(), 2)
        gene_result, var_result = get_info(replacement_responder)
        print("\nGenes found in multiple replacement/responder:")
        print_result(gene_result, dict(), 2)

    elif choice == "2":
        print("\nComparing Responders and non responders")
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"]['patient'])))

        gene_result_responder, var_result_responder = get_info(responder_list)
        gene_result_relapse, var_result_relapse = get_info(relapse_list)

        print("\nGenes found in multiple responders and no nonresponders:")
        print_result(gene_result_responder, gene_result_relapse, 2)
        print("\nSpecific variants found in multiple responders and no nonresponders:")
        print_result(var_result_responder, var_result_relapse, 3)
        print("\nGenes found in multiple nonresponders and no responders:")
        print_result(gene_result_relapse, gene_result_responder, 2)
        print("\nSpecific variants found in multiple nonresponders and no responders:")
        print_result(var_result_relapse, var_result_responder, 3)

    elif choice == "3":
        print("\nComparing evolution and no evolution")
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "evolution"]['patient'])))
        no_evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"]['patient'])))

        gene_result_evolution, var_result_evolution = get_info(evolution_list)
        gene_result_no_evolution, var_result_no_evolution = get_info(no_evolution_list)

        print("\nGenes found in multiple patients with evolution:")
        print_result(gene_result_evolution, gene_result_no_evolution, 2)
        print("\nSpecific variants found in multiple patients with evolution:")
        print_result(var_result_evolution, var_result_no_evolution, 3)
        print("\nGenes found in multiple patients without evolution:")
        print_result(gene_result_no_evolution, gene_result_evolution, 2)
        print("\nSpecific variants found in multiple patients without evolution:")
        print_result(var_result_no_evolution, var_result_evolution, 3)

    elif choice == "4":
        print("Comparing Responders and non responders within ibrutinib patients only")
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"].loc[df['treatment'] == "ibrutinib"]['patient'])))

        gene_result_responder, var_result_responder = get_info(responder_list)
        gene_result_relapse, var_result_relapse = get_info(relapse_list)

        print("\nGenes found in multiple responders:")
        print_result(gene_result_responder, gene_result_relapse, 2)
        print("\nSpecific variants found in multiple responders:")
        print_result(var_result_responder, var_result_relapse, 3)
        print("\nGenes found in multiple nonresponders:")
        print_result(gene_result_relapse, gene_result_responder, 2)
        print("\nSpecific variants found in multiple nonresponders:")
        print_result(var_result_relapse, var_result_responder, 3)

    elif choice == "5":
        print("\nComparing evolution and no evolution within ibrutinib patients only")
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "evolution"].loc[df['treatment'] == "ibrutinib"]['patient'])))
        no_evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"].loc[df['treatment'] == "ibrutinib"]['patient'])))

        gene_result_evolution, var_result_evolution = get_info(evolution_list)
        gene_result_no_evolution, var_result_no_evolution = get_info(no_evolution_list)

        print("\nGenes found in multiple patients with evolution:")
        print_result(gene_result_evolution, gene_result_no_evolution, 2)
        print("\nSpecific variants found in multiple patients with evolution:")
        print_result(var_result_evolution, var_result_no_evolution, 3)
        print("\nGenes found in multiple patients without evolution:")
        print_result(gene_result_no_evolution, gene_result_evolution, 2)
        print("\nSpecific variants found in multiple patients without evolution:")
        print_result(var_result_no_evolution, var_result_evolution, 3)

    elif choice == "6":
        print("\nComparing Responders and non responders within acalabrutinib patients only")
        responder_list = list(map(str,list(df.loc[df['response'] == "responder"].loc[df['treatment'] == "acalabrutinib"]['patient'])))
        relapse_list = list(map(str,list(df.loc[df['response'] == "relapse"].loc[df['treatment'] == "acalabrutinib"]['patient'])))

        gene_result_responder, var_result_responder = get_info(responder_list)
        gene_result_relapse, var_result_relapse = get_info(relapse_list)

        print("\nGenes found in multiple responders:")
        print_result(gene_result_responder, gene_result_relapse, 2)
        print("\nSpecific variants found in multiple responders:")
        print_result(var_result_responder, var_result_relapse, 3)
        print("\nGenes found in multiple nonresponders:")
        print_result(gene_result_relapse, gene_result_responder, 2)
        print("\nSpecific variants found in multiple nonresponders:")
        print_result(var_result_relapse, var_result_responder, 3)

    elif choice == "7":
        print("\nComparing evolution and no evolution within acalabrutinib patients only")
        evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "evolution"].loc[df['treatment'] == "acalabrutinib"]['patient'])))
        no_evolution_list = list(map(str,list(df.loc[df['evolution_status'] == "no_evolution"].loc[df['treatment'] == "acalabrutinib"]['patient'])))

        gene_result_evolution, var_result_evolution = get_info(evolution_list)
        gene_result_no_evolution, var_result_no_evolution = get_info(no_evolution_list)

        print("\nGenes found in multiple patients with evolution:")
        print_result(gene_result_evolution, gene_result_no_evolution, 2)
        print("\nSpecific variants found in multiple patients with evolution:")
        print_result(var_result_evolution, var_result_no_evolution, 3)
        print("\nGenes found in multiple patients without evolution:")
        print_result(gene_result_no_evolution, gene_result_evolution, 2)
        print("\nSpecific variants found in multiple patients without evolution:")
        print_result(var_result_no_evolution, var_result_evolution, 3)

    else:
        print("That was not an option, try again")
        choice = "enter"