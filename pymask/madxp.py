import itertools

import numpy as np
import pandas as pd

from cpymad.madx import Madx
import cpymad

def knob_df(my_knob,my_df):
    '''
    Filter the pandas DF, 'my_df', returning only the rows that depend on the selected knob, 'my_knob'.

    Args:
        my_knob: the name of the knob to filter.
        my_df: a pandas DF (it assumes that DF has a column called "knobs").
    Returns:
        The filter pandas DF showing the rows that depend on the selected knob, 'my_knob'.

    See madxp/examples/variablesExamples/000_run.py
    '''
    return my_df[my_df['knobs'].apply(lambda x: my_knob in x)]

def _knobs_from_parameters(parameters, indep_df, dep_df):
    '''
    Extract the list of knobs from a list of parameters.

    Args:
        parameters: list of parameters
        indep_df: independent variable DF
        dep_df: dependent variable DF
    Returns:
        The list of knobs corresponding to the list of parameters.

    See madxp/examples/variablesExamples/000_run.py
    '''
    my_knobs=[]
    for i in parameters:
        if i in indep_df.index:
            if not indep_df.loc[i]['constant']:
                my_knobs.append([i])
        else:
            try:
                my_knobs.append(dep_df.loc[i]['knobs'])
            except:
                print(f'Variable {i} not defined! Cosidered as a knob.')
                my_knobs.append([i])
    return list(itertools.chain.from_iterable(my_knobs))

def _extract_parameters(my_string):
    '''
    Extract all the parameters of a MAD-X expression.
    Args:
        my_string: The string of the MAD-X expression to parse.

    Returns:
        The list of the parameters present in the MAD-X expression.
    '''
    if (type(my_string)=='NoneType' or my_string==None
            or my_string=='None' or my_string=='[None]' or 'table(' in my_string):
        return []
    else:
        for i in [
        '*','->','-','/','+','^','(',')','[',']',',','\'','None']:
            my_string=my_string.replace(i,' ')
        my_list=my_string.split(' ')
        my_list=list(np.unique(my_list))
        if '' in my_list:
            my_list.remove('')
        if type(my_list)=='NoneType':
            my_list=[]
        for i in my_list.copy():
            if i.isdigit() or i[0].isdigit() or i[0]=='.':
                my_list.remove(i)
        my_list=list(set(my_list)-
        set([
            'sqrt',
            'log',
            'log10',
            'exp',
            'sin',
            'cos',
            'tan',
            'asin',
            'acos',
            'atan',
            'sinh',
            'cosh',
            'tanh',
            'sinc',
            'abs',
            'erf',
            'erfc',
            'floor',
            'ceil',
            'round',
            'frac',
            'ranf',
            'gauss',
            'tgauss']))
        return my_list


class Madxp(Madx):
    pass

    def set_variables_from_dict(self, params):
        for nn in params.keys():
            self.input(f'{nn}={params[nn]};')

    def get_sequences_df(self):
        '''
        Extract the pandas DF with the list of the sequences defined in the MAD-X handle.

        Returns:
            The pandas DF of the sequences. It can be and empty DF.

        See madxp/examples/variablesExamples/000_run.py
        '''
        sequences = self.sequence
        seq_dict = {}
        for ii in sequences:
            seq_dict[ii] = {}
            if sequences[ii].has_beam:
                seq_dict[ii]['beam'] = True
            else:
                seq_dict[ii]['beam'] = False
            if sequences[ii].is_expanded:
                seq_dict[ii]['expanded'] = True
            else:
                seq_dict[ii]['expanded'] = False
        return pd.DataFrame(seq_dict).transpose()

    def get_beams_df(self):
        '''
        Extract the pandas DF with the beams associated to the sequences defined in the MAD-X handle.

        Returns:
            The pandas DF of the beams. It can be and empty DF.
        See madxp/examples/variablesExamples/000_run.py
        '''
        df_list = []
        sequences = self.sequence
        for ii in sequences:
            try:
                df_list.append(pd.DataFrame([dict(sequences[ii].beam)], index=[ii]))
            except:
                print(f'The sequence {ii} has no beam attached.')
        if len(df_list) > 0:
            return pd.concat(df_list)
        else:
            return pd.DataFrame()

    def get_variables_dicts(self, expressions_as_str=True):
        variables_df = self.get_variables_dataframes()
        outp = {
            'constants': variables_df['constants'].to_dict()['value'],
            'independent_variables':
                variables_df['independent_variables'].to_dict()['value'],
            'dependent_variables_expr':
                variables_df['dependent_variables'].to_dict()['expression'],
            'dependent_variables_val':
                variables_df['dependent_variables'].to_dict()['value'],
            }

        outp['all_variables_val'] = {kk:outp['constants'][kk] for
                kk in outp['constants'].keys()}
        outp['all_variables_val'].update(outp['independent_variables'])
        outp['all_variables_val'].update(outp['dependent_variables_val'])


        return outp

    def get_variables_dataframes(self, expressions_as_str=True):
        '''
        Extract the dictionary of the variables and constant pandas DF of the -X global workspace.

        Returns:
            The a dictionary containing:
            - constants_df: the pandas DF with the constants
            - independent_variables: the pandas DF with the independent variables
            - dependent_variables: the pandas DF with the dependent variables
        All the three DFs have a columns 'value' with the numerical values of the costants/variables.
        The dependent_variable_df, in addition to 'value' has the following columns:
            - 'expression': the string corrensponding to the MAD-X expression
            - 'parameters': the list of parameters used in the expression
            - 'knobs': the list of the independent variables that control
              the dependent variables. Note tha the parameters can be constants and/or dependent variables,
              whereas the 'knobs' are only independent variables.
        '''
        my_dict={}
        aux=self._independent_variables_df()
        import numpy as np
        independent_variables_df=aux[np.logical_not(aux['constant'])].copy()
        del independent_variables_df['constant']
        constant_df=aux[aux['constant']].copy()
        del constant_df['constant']
        my_dict['constants']=constant_df
        my_dict['independent_variables']=independent_variables_df
        my_dict['dependent_variables']=self._dependent_variables_df()

        if expressions_as_str:
            my_dict['dependent_variables']['expression'] = (
                     my_dict['dependent_variables']['expression'].apply(
                         str))
        return my_dict

    def _dependent_variables_df(self):
        '''
        Extract the pandas DF with the dependent variables of the MAD-X handle.

        Returns:
            The pandas DF of the dependent variables. The columns of the DF correspond to the
            - the numerical value of the dependent variable (value)
            - the string corrensponding to the MAD-X expression (expression)
            - the list of parameters used in the expression (parameters)
            - the list of the fundamental independent variables.
              These are independent variables that control numerical values of the variable (knobs).

        See madxp/examples/variablesExamples/000_run.py
        '''
        my_dict={}
        for i in list(self.globals):
            aux=_extract_parameters(str(self._libmadx.get_var(i)))
            if aux!=[]:
                my_dict[i]={}
                my_dict[i]['parameters']=list(np.unique(aux))

        myhash=hash(str(my_dict))
        while True:
            for i in my_dict:
                aux=[]
                for j in my_dict[i]['parameters']:
                    try:
                        aux=aux+my_dict[j]['parameters']
                    except:
                        aux=aux+[j]
                my_dict[i]['knobs']=list(np.unique(aux))
            if myhash==hash(str(my_dict)):
                break
            else:
                myhash=hash(str(my_dict))

        for i in my_dict:
            for j in my_dict[i]['knobs'].copy():
                if self._libmadx.get_var_type(j)==0:
                    my_dict[i]['knobs'].remove(j)
            my_dict[i]['expression']=self._libmadx.get_var(i)
            my_dict[i]['value']=self.globals[i]

        if len(my_dict)>0:
            return pd.DataFrame(my_dict).transpose()[['value','expression','parameters','knobs']].sort_index()
        else:
            return pd.DataFrame()

    def _independent_variables_df(self):
        '''
        Extract the pandas DF with the independent variables of the MAD-X handle.

        Returns:
            The pandas DF of the independent variables. The columns of the DF correspond to the
            - the numerical value of the independent variable (value)
            - a boolean value to know it the variable is constant or not (constant)

        See madxp/examples/variablesExamples/000_run.py
        '''

        dep_df=self._dependent_variables_df()
        #if len(dep_df)>0:
        #    aux=list(dep_df['knobs'].values)
        #    aux=list(itertools.chain.from_iterable(aux))
        #fundamentalSet=set(np.unique(aux))
        independent_variable_set=set(self.globals)-set(dep_df.index)
        my_dict={}
        for i in independent_variable_set:
            my_dict[i]={}
            if self._libmadx.get_var_type(i)==0:
                my_dict[i]['constant']=True
                # my_dict[i]['knob']=False
            else:
                my_dict[i]['constant']=False
                # if i in fundamentalSet:
                #    my_dict[i]['knob']=True
                # else:
                #    my_dict[i]['knob']=False
            my_dict[i]['value']=self.globals[i]

        return pd.DataFrame(my_dict).transpose()[['value','constant']].sort_index()


    def get_sequence_df(self,sequenceName):
        '''
        Extract a pandas DF of the list of the elements and all their attributes for a given sequence.

        Args:
            sequenceName: the sequence name
        Returns:
            The list of knobs corresponding to the list of parameters.

        See madxp/examples/variablesExamples/000_run.py
        '''
        my_list=[]
        sequences=self.sequence
        my_sequence=sequences[sequenceName]
        indep_df=self._independent_variables_df()
        dep_df=self._dependent_variables_df()

        for my_index, _ in enumerate(my_sequence.elements):
            aux=self._libmadx.get_element(sequenceName,my_index)
            my_dict={}
            my_dict['parameters']=[]
            for i in aux:
                my_dict[i]=aux[i]
            del my_dict['data']

            for i in aux['data']:
                my_dict[i]=str(aux['data'][i])
                if isinstance(aux['data'][i], cpymad.types.Parameter):
                    my_dict[i+' value']=aux['data'][i].value
                    my_dict['parameters']+=_extract_parameters(str(aux['data'][i].expr))
            my_dict['parameters']=np.unique(my_dict['parameters'])
            my_list.append(my_dict)
        my_df=pd.DataFrame(my_list)
        my_df=my_df.set_index('name')
        my_df.index.name=''
        my_df['knobs']=my_df['parameters'].apply(lambda x: _knobs_from_parameters(x,indep_df,dep_df))
        first_columns=['position','parent','base_type','length','parameters','knobs']
        last_columns=list(set(my_df.columns)-set(first_columns))
        last_columns.sort()
        return my_df[first_columns+last_columns]

    def get_twiss_df(self, table_name):
        '''
        Extract the pandas DF of a MAD-X table.

        Args:
            table_name: Name of the table

        Returns:
            The pandas DF of a MAD-X table.

        See madxp/examples/variablesExamples/000_run.py
        '''
        table = self.table[table_name]
        my_df=pd.DataFrame(dict(table))
        my_df=my_df.set_index('name', drop = False)
        my_df.index.name=''
        return my_df

    def get_summ_df(self, table_name):
        '''
        Extract the pandas DF of a MAD-X summary table.

        Args:
            table_name: Name of the table
        Returns:
            The pandas DF of a MAD-X table.

        See madxp/examples/variablesExamples/000_run.py
        '''
        table = self.table[table_name]
        my_df=pd.DataFrame(dict(table))
        my_df.index=[table._name]
        return my_df

