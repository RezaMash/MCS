import sys
import ast
sys.path.insert(0, '/home/rmiraska/Workspace/MetabolicNetworks/MongooseGUI3')
import ModelParsing
import ModelProcessing

path_to_models = '/home/rmiraska/Workspace/MetabolicNetworks/models/'


def fraction_to_floating_point(vector):
    return [float(i) for i in vector]


def round_a_vector(vector, ndigits):
    return [round(i, ndigits) for i in vector]


def readable(float_number):
    answer = str(float_number)
    if answer.endswith('.0'):
        answer = answer[:-2]
    return answer.rjust(8, ' ')


def sbml_to_txt(model_name):
    network = ModelParsing.parseSBML(path_to_models+model_name+'/sbml.xml')
    full_matrix = network.fullMatrix

    irreversible = network.findIrreversibleReactions()
    rev = [0 if i+1 in irreversible else 1 for i in range(len(network.reactions))]
    integralized_full_matrix = [ModelProcessing.integralize(row) for row in full_matrix]
    floated_full_matrix = [fraction_to_floating_point(row) for row in full_matrix]
    rounded_full_matrix = [round_a_vector(row, 2) for row in floated_full_matrix]
    open(path_to_models+model_name+'/full_matrix_python_format.txt', 'w').write(str(full_matrix))
    open(path_to_models + model_name + '/reversibility_python_format.txt', 'w').write(str(rev))
    open(path_to_models + model_name + '/full_matrix_integralized_python_format.txt', 'w')\
        .write(str(integralized_full_matrix))
    open(path_to_models + model_name + '/full_matrix_floated_python_format.txt', 'w')\
        .write(str(floated_full_matrix))
    file = open(path_to_models + model_name + '/full_matrix_human_readable.txt', 'w')
    file.write(str(len(full_matrix[0]))+'\n'+str(len(full_matrix))+'\n')
    for row in rounded_full_matrix:
        for i in row:
            file.write(readable(i))
        file.write('\n')

    network.reduceNetwork()
    reduced_matrix = network.reducedMatrix
    reduced_rev = [1 if network.reactionSubsets[i].reversible else 0 for i in range(len(network.reducedMatrix[0]))]

    integralized_reduced_matrix = [ModelProcessing.integralize(row) for row in reduced_matrix]
    floated_reduced_matrix = [fraction_to_floating_point(row) for row in reduced_matrix]
    rounded_reduced_matrix = [round_a_vector(row, 2) for row in floated_reduced_matrix]
    open(path_to_models + model_name + '/reduced_matrix_python_format.txt', 'w').write(str(reduced_matrix))
    open(path_to_models + model_name + '/reduced_reversibility_python_format.txt', 'w').write(str(reduced_rev))
    open(path_to_models + model_name + '/reduced_matrix_integralized_python_format.txt', 'w') \
        .write(str(integralized_reduced_matrix))
    open(path_to_models + model_name + '/reduced_matrix_floated_python_format.txt', 'w') \
        .write(str(floated_reduced_matrix))
    file = open(path_to_models + model_name + '/reduced_matrix_human_readable.txt', 'w')
    file.write(str(len(reduced_matrix[0])) + '\n' + str(len(reduced_matrix)) + '\n')
    for row in rounded_reduced_matrix:
        for i in row:
            file.write(readable(i))
        file.write('\n')



def txt_to_dlm(model_name):
    matrix = open(path_to_models + model_name +'/full_matrix_integralized_python_format.txt', 'r').read()
    matrix = matrix.replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '; ')
    open(path_to_models + model_name + '/full_matrix_integralized_dlm_matlab_format.txt', 'w').write(matrix)

    matrix = open(path_to_models + model_name + '/full_matrix_floated_python_format.txt', 'r').read()
    matrix = matrix.replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '\n ')
    open(path_to_models + model_name + '/full_matrix_floated_dlm_matlab_format.txt', 'w').write(matrix)

    rev = open(path_to_models + model_name + '/reversibility_python_format.txt', 'r').read()
    rev = rev.replace('[', '')
    rev = rev.replace(']', '')
    open(path_to_models + model_name + '/reversibility_dlm_matlab_format.txt', 'w').write(rev)

    matrix = open(path_to_models + model_name + '/reduced_matrix_integralized_python_format.txt', 'r').read()
    matrix = matrix.replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '; ')
    open(path_to_models + model_name + '/reduce_matrix_integralized_dlm_matlab_format.txt', 'w').write(matrix)

    matrix = open(path_to_models + model_name + '/reduced_matrix_floated_python_format.txt', 'r').read()
    matrix = matrix.replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '\n ')
    open(path_to_models + model_name + '/reduced_matrix_floated_dlm_matlab_format.txt', 'w').write(matrix)

    rev = open(path_to_models + model_name + '/reduced_reversibility_python_format.txt', 'r').read()
    rev = rev.replace('[', '')
    rev = rev.replace(']', '')
    open(path_to_models + model_name + '/reduced_reversibility_dlm_matlab_format.txt', 'w').write(rev)


def line_to_vector(line):
    for i in range(10):
        line = line.replace('  ', ' ')

    if line.startswith(' '):
        line = line[1:]

    if line.endswith(' '):
        line = line[:-1]

    # for i in line.split(' '):
    #     print('!'+i)

    return [float(element) for element in line.split(' ')]


# def hr_to_txt(model_name):
#     file = open(path_to_models + model_name + '/full_matrix_human_readable.txt', 'r')
#     full_matrix = [line_to_vector(line) for line in iter(lambda: file.readline(), '')]
#     open(path_to_models + model_name + '/full_matrix_floated_python_format.txt', 'w') \
#         .write(str(full_matrix))
#     rev = open(path_to_models + model_name + '/reversibility_dlm_matlab_format.txt', 'r').read()
#     if rev.endswith('\n'):
#         rev = rev[:-1]
#     rev = '[' + rev + ']'
#     open(path_to_models + model_name + '/reversibility_python_format.txt', 'w') \
#         .write(rev)


def sbml_to_null(model_name):
    network = ModelParsing.parseSBML(path_to_models + model_name + '/sbml.xml')
    null_basis = ModelProcessing.NullspaceBasis(network.fullMatrix)
    null_basis = ModelProcessing.transpose(null_basis)

    integralized_null_basis = [ModelProcessing.integralize(row) for row in null_basis]
    floated_null_basis = [fraction_to_floating_point(row) for row in null_basis]
    rounded_null_basis = [round_a_vector(row, 2) for row in floated_null_basis]
    open(path_to_models + model_name + '/null_basis_python_format.txt', 'w').write(str(null_basis))
    open(path_to_models + model_name + '/null_basis_integralized_python_format.txt', 'w') \
        .write(str(integralized_null_basis))
    open(path_to_models + model_name + '/null_basis_floated_python_format.txt', 'w') \
        .write(str(floated_null_basis))

    matrix = str(floated_null_basis).replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '\n ')
    open(path_to_models + model_name + '/null_basis_floated_dlm_matlab_format.txt', 'w').write(matrix)

    file = open(path_to_models + model_name + '/null_basis_human_readable.txt', 'w')
    for row in rounded_null_basis:
        for i in row:
            file.write(readable(i))
        file.write('\n')

    network.reduceNetwork()
    reduced_null_basis = ModelProcessing.NullspaceBasis(network.reducedMatrix)
    reduced_null_basis = ModelProcessing.transpose(reduced_null_basis)

    integralized_reduced_null_basis = [ModelProcessing.integralize(row) for row in reduced_null_basis]
    floated_reduced_null_basis = [fraction_to_floating_point(row) for row in reduced_null_basis]
    rounded_reduced_null_basis = [round_a_vector(row, 2) for row in floated_reduced_null_basis]
    open(path_to_models + model_name + '/reduced_null_basis_python_format.txt', 'w').write(str(reduced_null_basis))
    open(path_to_models + model_name + '/reduced_null_basis_integralized_python_format.txt', 'w') \
        .write(str(integralized_reduced_null_basis))
    open(path_to_models + model_name + '/reduced_null_basis_floated_python_format.txt', 'w') \
        .write(str(floated_reduced_null_basis))

    matrix = str(floated_reduced_null_basis).replace('[[', '')
    matrix = matrix.replace(']]', '')
    matrix = matrix.replace('], [', '\n ')
    open(path_to_models + model_name + '/reduced_null_basis_floated_dlm_matlab_format.txt', 'w').write(matrix)

    file = open(path_to_models + model_name + '/reduced_null_basis_human_readable.txt', 'w')
    for row in rounded_reduced_null_basis:
        for i in row:
            file.write(readable(i))
        file.write('\n')


for model in ['e_coli_core', 'BIOMD0000000001', 'BIOMD0000000002', 'BIOMD0000000003', 'BIOMD0000000004',
              'BIOMD0000000005', 'BIOMD000000034', 'BIOMD0000000610', 'BIOMD0000000674']:
    sbml_to_txt(model)
    txt_to_dlm(model)
    sbml_to_null(model)

# hr_to_txt('simple_01')