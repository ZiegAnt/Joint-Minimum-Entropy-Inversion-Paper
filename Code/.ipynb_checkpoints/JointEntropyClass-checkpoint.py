# %% 
# Modified Version that allows differing start and reference models


# %% Imports
import  pygimli                 as     pg
import  numpy                   as     np
import  matplotlib.pyplot       as     plt
import  sys
# %% JointEntropyInversion class
class JointEntropyInversion(object):
    
    def __init__(self, mgr_list=None, data_list=None, mesh=None, beta=1e-10, order=1, smooth_factor=1):
        '''
        Initiating Joint Entropy Inversion class after Zhdanov et al (2022)

        Parameters
        ----------
        mgr_list : list, optional
            List of method managers. Can also be set via setManagers. The default is None.
        data_list : list, optional
            List of data. Can also be set via setData. The default is None.
        mesh : pg.mesh, optional
            Inversion mesh that is building the paraDomain for every method. 
            Can also be set via setMesh. The default is None.
        beta : float, optional
            Numerical stabilizer. SHould be much smaller than model parameters. 
            The default is 1e-10.
        order : float or int, optional
            order of entropy representation. Should not be higher than 2. 
            The default is 1.
        smooth_factor : list or float or int, optional
            smoothing factors that determines the smoothing strength w.r.t. the
            mean modelweight value for each method. Default is 1. 
            If float or int is given then smoothing factor is same for all methods

        Returns
        -------
        None.

        '''        
        # Methods and data
        self.mgrs      = mgr_list
        self.data_list = data_list
        self.names     = [f'Method_{i}' for i in range(len(self.mgrs))]
        
        # Mesh/Grid
        self.mesh      = mesh
        self.mesh_list = [self.mesh for m in self.mgrs]
        
        # Parameters
        self.q = order         # order of entropy expression
        self.b = beta          # numerical stabilizer
        self.method_weights = [1 for i in range(len(self.mgrs))] # method specific weighting 
        
        if type(smooth_factor)==list:
            self.a = smooth_factor         # factor that determines smoothing w.r.t. to median of mWeights; a=0 indicates pure damping
        else:
            self.a = [smooth_factor for i in range(len(self.mgrs))]
        
        # Temporary helper variables
        self.Q_JME    = 0
        self.Q_JMEG   = 0
        self.Grad_sum = 0
        self.Diff_sum = 0
        
        # Inversion Keyword dicts
        self.KWInv = [None for m in self.mgrs]
        
        # Reference models
        self.refmod = [None for m in self.mgrs]
        
        # Saving variables
        self.weightsHistory    = 'There was no inversion run so there is no history yet'
        self.modelsHistory     = 'There was no inversion run so there is no history yet'
        self.responseHistory     = 'There was no inversion run so there is no history yet'
        self.chi2History       = 'There was no inversion run so there is no history yet'
        self.phiHistory        = 'There was no inversion run so there is no history yet'
        self.phiMHistory       = 'There was no inversion run so there is no history yet'
        self.phiDHistory       = 'There was no inversion run so there is no history yet'
        self.rmsHistory        = 'There was no inversion run so there is no history yet'
        self.jointChi2History  = 'There was no inversion run so there is no history yet'
        self.SHistory          = 'There was no inversion run so there is no history yet'
        self.log_JME_History   = 'There was no inversion run so there is no history yet'
        self.log_JMEG_History  = 'There was no inversion run so there is no history yet'
        
        
    def setParameters(self, q=1, b=1e-10, a=1):
        '''
        Set Joint inversion parameters.
        Parameters
        ----------
        q : float or int
            order for min. entropy constraints. The default is 1.
        b : float
            numerical stabilization factor. Should be much smaller than model 
            parameter values. The default is 1e-10.
        a : list or float or int
            smoothing factors that determines the smoothing strength w.r.t. the
            mean modelweight value for each method. Default is 1. 
            If float or int is given then smoothing factor is same for all methods
        '''
        self.q = q
        self.b = b
        if type(a)==list:
            self.a = a
        else:
            self.a = [a for i in range(len(self.mgrs))]

    def setMesh(self, mesh, for_all=True):
        '''
        Set mesh for inversion constraint calculations
        Parameters
        ----------
        mesh : pg.mesh
        for_all : boolean
            If true then the mesh is set as inversion mesh for all methods. 
            Note that ERT inversions should have an additional boundary region.
            The default is True.
        '''
        self.mesh = mesh
        if for_all:
            for i in range(len(self.mgrs)):
                self.setMethodMesh(mesh,i)

    def setMethodMesh(self, mesh, method_id):
        """
        Set Inversion mesh or grid for specific method.
        Thought for ERT grids that hold an additional boundary region for inversion
        """
        if self.mesh_list ==[]:
            self.mesh_list = [self.mesh for m in self.mgrs]
            
        self.mesh_list[method_id] = mesh    
        
    def setMethodWeights(self, weightList):
        '''
        Set list of weighting factors that weights the method w.r.t. eachother.
        '''
        if len(weightList)==len(self.mgrs):
            self.method_weights = weightList
        else:
            raise pg.critical("Number of given method weight factors does not match methods in self.mgrs!!! Check again!!!")
        

    def setData(self, data_list):
        """
        Set list of data.
        """
        if len(self.mgrs) != len(data_list):
            self.data_list = data_list
            for i in range(len(self.mgrs)):
                self.mgrs[i].setData(self.data_list[i])
        else:    
            raise pg.critical("Length of data list not matching number of methods")
            
    def setManagers(self, mgr_list, names=None):
        """
        Set list of Managers. 
        Optional add names for each method by including corresponding list.
        """
        self.mgrs = mgr_list
        if names!=None:
            self.setNames(names)
    
    def setNames(self, names):
        """
        Set list of Method names.
        """
        self.names = names

    def setKWInv(self, dict_list):
        '''
        Set inversion keyword dictionaries.
        Possible entries could be:
            lam
            cType (should be 0 or 10)
            lambdaFactor
            verbose
            zWeight
            startModel
            vTop (for SRT)
            vBottom (for SRT)
            etc. 
        see pygimli documentation for more inversion parameters

        Parameters
        ----------
        dict_list : list
            List containing one kwarg dictionary for each method.
        '''
        self.KWInv = dict_list
        
    def setReferenceModel(self,refmod_list):
        """
        Set list of Reference Models (list of pg.Vectors).
        """
        self.refmod = refmod_list

    def normalize(self, m, m0):
        '''
        Transform starting model and current model to values between 0 and 1 by
        doing linear normalization of logarithmic model. This ensures comparability of different models

        Parameters
        ----------
        m : np.array or pg.Vector
            current model.
        m0 : np.array or pg.Vector
            starting model.

        Returns
        -------
        m_new, m0_new: normalized starting and current model for weight calculations

        '''
        mmin = min(np.log(m))
        mmax = max(np.log(m))
        m_new = (np.log(m)-mmin)/(mmax-mmin+self.b)
        m0_new = (np.log(m0)-mmin)/(mmax-mmin+self.b)
        return m_new, m0_new
        

    def calcQ_JME(self):
        '''
        Calculate Volume Integral of all models.
        Used for normaization in min. Entropy Model Weighting function.
        '''
        Q = []
        for mgr_i, mgr in enumerate(self.mgrs):
            m0 = self.refmod[mgr_i]
            m = mgr.model
            
            [m, m0] = self.normalize(m, m0)
            
            Qi = np.sum(self.method_weights[mgr_i]*abs(m-m0)**self.q + self.b)
            Q.append(Qi)   
        return sum(Q)
    
    def calcDiff_sum_vec(self):
        '''
        Calculate vector (dimension of model space) containing sum of 
        differences to reference models for all methods. 
        (Denominator in log function for model weights)
        '''        
        diff = np.zeros((self.mesh.cellCount(),len(self.mgrs)))
        
        for i in range(len(self.mgrs)):
            m0 = self.refmod[i]
            m = self.mgrs[i].model
            
            [m, m0] = self.normalize(m, m0)
            
            diff[:,i] = self.method_weights[i]*abs(m-m0)**self.q
        return np.sum(diff, axis=1)
    
    def calcQ_JMEG(self):
        '''
        Calculate Volume Integral of all models.
        Used for normaization in min. Entropy gradient Model Weighting function.
        '''
        Q = []
        for mgr_i, mgr in enumerate(self.mgrs):
            m0 = self.refmod[mgr_i]
            m = mgr.model
            [m, m0] = self.normalize(m, m0)   
            grad_m = pg.solver.grad(self.mesh, m)
            
            Qi = np.sum(np.linalg.norm(grad_m,axis=1)**self.q + self.b)
            Q.append(Qi)
        return sum(Q)
    
    def calcGrad_sum_vec(self):
        '''
        Calculate vector (dimension of model space) containing sum of 
        gradients for all methods. 
        (Denominator in log function for model weights)
        '''   
        grad = np.zeros((self.mesh.cellCount(),len(self.mgrs)))
        
        for i in range(len(self.mgrs)):
            m0 = self.refmod[i]
            m = self.mgrs[i].model
            [m, m0] = self.normalize(m, m0)
            grad_m = pg.solver.grad(self.mesh, m)
            
            grad[:,i] = self.method_weights[i]*np.linalg.norm(grad_m,axis=1)**self.q
        return np.sum(grad, axis=1)
    
        
    def setMWeightMinEntropy(self, iteration, inv):
        '''
        Calculate Model weights according to Zhdanov 2022 with minimum entropy 
        constraints and set them to cWeights in current inversion instance.
        Note: the smoothing constraints are set to a value np.mean(mWeights)*smooth_factor
        Parameters
        ----------
        iteration : Int
            Iteration number.
        inv : pygimli.frameworks.inversion.Inversion
            Inversion framework.
        '''
        if iteration>0:
            mgr_i = np.argwhere([inv == mgr.inv for mgr in self.mgrs])[0][0] # Index of manager in self.mgrs
            
            m0 = self.refmod[mgr_i]
            m = inv.model
            [m, m0] = self.normalize(m, m0)
            
            mesh = self.mesh
            Q = self.Q_JME
        
            we = np.sqrt( 1/Q * 1/(self.b+abs(m-m0)**2) * (self.b+self.method_weights[mgr_i]*abs(m-m0)**self.q) * np.log(Q/(self.b+self.Diff_sum)))
            
            Cweights = np.ones(len(inv.inv.cWeight()))*np.mean(we)*self.a[mgr_i] # re-set constrain weights
            Cweights[-mesh.cellCount():] = we                                    # re-set constrain weights corresponding to model weights
            Cweights = pg.Vector(Cweights)
            inv.inv.setCWeight(Cweights)                                         # set weights in inversion instance


    def setMWeightMinEntropyGradient(self, iteration, inv):
        '''
        Calculate Model weights according to Zhdanov 2022 with minimum entropy 
        gradient constraints and set them to cWeights in current inversion instance.
        Note: the smoothing constraints are set to a value np.mean(mWeights)*smooth_factor

        Parameters
        ----------
        iteration : Int
            Iteration number.
        inv : pygimli.frameworks.inversion.Inversion
            Inversion framework.
        '''
        if iteration>0:
            mgr_i = np.argwhere([inv == mgr.inv for mgr in self.mgrs])[0][0] # Index of manager in self.mgr
            
            m0 = self.refmod[mgr_i]
            m = inv.model
            [m, m0] = self.normalize(m, m0)
            
            mesh = self.mesh
            grad_m = pg.solver.grad(self.mesh, m)     
            Q = self.Q_JMEG
            
            we = np.sqrt( 1/Q * 1/(self.b+abs(m-m0)**2) * (self.b+self.method_weights[mgr_i]*np.linalg.norm(grad_m,axis=1)**self.q) * np.log(Q/(self.b+self.Grad_sum)))
            
            Cweights = np.ones(len(inv.inv.cWeight()))*np.mean(we)*self.a[mgr_i] # re-set constrain weights 
            Cweights[-mesh.cellCount():] = we                                    # re-set constrain weights corresponding to model weights 
            Cweights = pg.Vector(Cweights)
            inv.inv.setCWeight(Cweights)                                         # set weights in inversion instance


    def runInversion(self, method ,maxIter=10, KWInv=None, breakup_criterion='joint', chi_limit=2):
        
        '''
        Run minimum entropy constraint joint inversion.
        Note: the reference model for the damping is always set to be the starting model. 

        Parameters
        ----------
        method : str
            'ME' for min. entropy constraint, 'MEG' for min. entropy gradient constraint.
        maxIter : int, optional
            maximum number of iterations. The default is 10.
        KWInv : list, optional
            List containing inversion kwarg dictionaries. The default is None.
        breakup_criterion : string, optional
            Type of breakup criterion. Either 'joint' for the joint Chi^2 criterion or 
            'all' for min. Chi^2 for all methods. The default is 'joint'.
        chi_limit : float, optional
            Break-ip criterion of inversion. The default is 2.
        '''        
        
        if KWInv!=None:
            self.setKWInv(KWInv)
            
        # Check if necessary keyword argumant dicts are defined for all methods
        if any(dic==None for dic in self.KWInv):
            raise pg.critical("Inversion KW dictionaries are not set for all methods!!! Set them via self.setKWInv")
            sys.exit()
            
        # Check if necessary reference models are defined for all methods
        if any(rm==None for rm in self.refmod):
            raise pg.critical("Reference models are not set for all methods!!! Set them via self.setReferenceModel")
            sys.exit()
        
        # Define saving variables
        self.weightsHistory   = [[] for mgr in self.mgrs]
        self.modelsHistory    = [[] for mgr in self.mgrs]
        self.responseHistory    = [[] for mgr in self.mgrs]
        self.chi2History      = [[] for mgr in self.mgrs]
        self.phiHistory       = [[] for mgr in self.mgrs]
        self.phiMHistory       = [[] for mgr in self.mgrs]
        self.phiDHistory       = [[] for mgr in self.mgrs]
        self.rmsHistory       = [[] for mgr in self.mgrs]
        self.jointChi2History = []
        self.SHistory         = []
        self.log_JME_History    = []
        self.log_JMEG_History   = []
        
        # Run one iteration for each method to set inversion instances up
        print('Starting Iteration 0...')
        
        for mgr_i,mgr in enumerate(self.mgrs):
            print(f'... {self.names[mgr_i]}')
            # Set reference models to inversion frame works
            mgr.inv.inv.setReferenceModel(self.refmod[mgr_i])
            # Invert first iteration to ensure that all helpers in pygimli are set properly
            mgr.invert(maxIter=1, mesh=self.mesh_list[mgr_i], isReference=False, cType=10, **self.KWInv[mgr_i])
            # This ensures that the inv.model of a traveltime manager is in slowness
            if type(mgr)==pg.physics.traveltime.TravelTimeManager:
                mgr.inv.model = 1/mgr.inv.model
        
        # Loop over iterations
        for i in range(self.mgrs[0].inv.inv.iter(), maxIter+1):
            joint_chi2 = []
            Si = []
            
            print(f'Starting iteration {i}...')
            
            # Define helper variables
            if method == 'ME':
                self.Q_JME = self.calcQ_JME()
                self.Diff_sum = self.calcDiff_sum_vec()
                self.log_JME_History.append(np.log(self.Q_JME/(self.b+self.Diff_sum)))
            elif method == 'MEG':
                self.Q_JMEG = self.calcQ_JMEG()
                self.Grad_sum = self.calcGrad_sum_vec()
                self.log_JMEG_History.append(np.log(self.Q_JMEG/(self.b+self.Grad_sum)))
                
            else:
                raise pg.critical("Method must be joint min. Entropy (ME) or joint min. Entropy Gradient (MEG)")
                sys.exit()
            
            # Set model weights
            for mgr_i,mgr in enumerate(self.mgrs):
                self.modelsHistory[mgr_i].append(mgr.inv.model)
                self.responseHistory[mgr_i].append(mgr.inv.response)
                print(f'Set J{method} weights for {self.names[mgr_i]} with q={self.q}, b={self.b} and a={self.a[mgr_i]}')
                if method == 'ME':
                    self.setMWeightMinEntropy(i, mgr.inv)
                elif method == 'MEG':
                    self.setMWeightMinEntropyGradient(i, mgr.inv)

                
            # run one iteration for each method    
            for mgr_i,mgr in enumerate(self.mgrs):    
                mgr.inv.inv.oneStep() 
                mgr.inv.model = mgr.inv.inv.model()
                
                # Save important values
                we = np.array(mgr.inv.inv.cWeight())[-self.mesh.cellCount():]
                self.weightsHistory[mgr_i].append(we)
                self.phiHistory[mgr_i].append(mgr.inv.phi())
                self.phiMHistory[mgr_i].append(mgr.inv.phiModel())
                self.phiDHistory[mgr_i].append(mgr.inv.phiData())
                self.rmsHistory[mgr_i].append(mgr.inv.relrms())
                self.chi2History[mgr_i].append(mgr.inv.chi2())
                joint_chi2.append(mgr.inv.chi2())
                
                # Entropy stabilizing functional
                m0 = pg.Vector(mgr.inv.startModel)
                m = mgr.inv.model
                [m, m0] = self.normalize(m, m0)
                
                Si.append(np.dot(we*(m-m0),we*(m-m0))) 
                
                print(f'#####     {self.names[mgr_i]} weighted misfit chi^2   = {mgr.inv.chi2():.2f}')

            self.jointChi2History.append(sum(joint_chi2))
            self.SHistory.append(sum(Si))
            
            print(f'#####     stabilizing functional    S = {sum(Si):.2f}')
            print(f'#####     joint weighted misfit chi^2 = {sum(joint_chi2):.2f}')
            print('#'*60)
            print('       ')
            
            if breakup_criterion=='joint' and sum(joint_chi2) <= chi_limit:
                print('~'*60)
                print(f'          joint chi^2 < {chi_limit:.2f}         ')
                print('~'*60)
                break
            
            elif breakup_criterion=='all' and all(mgr.inv.chi2() <= chi_limit for mgr in self.mgrs):
                print('~'*60)
                for mgr_i,mgr in enumerate(self.mgrs):
                    print(f'            {self.names[mgr_i]} weighted misfit chi^2   = {mgr.inv.chi2():.2f}')
                print('~'*60)
                break
            
        self.weightsHistory = np.array(self.weightsHistory)
        print('All weights saved under self.weights')
        
        for mgr_i,mgr in enumerate(self.mgrs):
            mgr.inv.chi2History = self.chi2History[mgr_i]
    
    def getModels(self):
        '''
        Returns
        -------
        m_est : list
            List containing the final model parameter vectors.
        '''
        m_est = []
        for mgr in self.mgrs:
            m_est.append(mgr.model)
        return m_est

    def plotFitHistory(self):
        
        n = len(self.mgrs) # number of methods
        it = len(self.chi2History[0]) # number of iterations that were done
        
        fig, ax = plt.subplots(1,n+1, figsize=(12,4))
        fig.tight_layout(pad=6)
        for i in range(n):

            ax[i].plot(np.arange(it)+1, self.chi2History[i], 'r-o')
            ax[n].plot(np.arange(it)+1, self.chi2History[i], '-o', label=self.names[i])
            ax[i].set_xlabel('Iteration')
            ax[i].set_ylabel('$\chi^2$', color='red')
            ax[i].tick_params(axis ='y', which='both', labelcolor='red', labelsize=8)
            ax[i].set_title(f'Fit for {self.names[i]}')
            ax[i].grid(which='minor',axis='y', color='red', linestyle=':', linewidth=0.4)
            ax[i].set_yscale('log')

            ax2 = ax[i].twinx()
            
            ax2.plot(np.arange(it)+1, self.rmsHistory[i], 'b-*')
            ax2.set_ylabel('$RMS-Error$ in $\%$', color='blue')
            ax2.tick_params(axis ='y', which='both', labelcolor='blue', labelsize=8)
            ax2.grid(which='minor',axis='y', color='blue', linestyle=':', linewidth=0.4)

        ax[n].set_xlabel('Iteration')
        ax[n].set_ylabel('$\chi^2$')
        ax[n].set_title('$\chi^2$ comparison')
        ax[n].set_yscale('log')
        ax[n].tick_params(axis ='y', which='both',labelsize=8)
        ax[n].grid(which='minor',axis='y',linestyle=':', linewidth=0.4)
        ax[n].legend(loc='upper right')
        
        
    def plotWeights(self, method, step=4):
        
        n = len(self.mgrs) # number of methods
        it = len(self.chi2History[0]) # number of iterations that were done
        nr = int(it/step) #number of screenshots
        
        fig, ax = plt.subplots(n+1,nr, figsize=(12,4))
        # fig.tight_layout(pad=6)
        
        if method=='ME':
            d = self.log_JME_History
            tit = 'ln(Q / (Diffsum+beta) '
        else:
            d = self.log_JMEG_History
            tit = 'ln(Q / (Gradsum+beta) '
        
        for j in range(nr):
            
            pg.show(self.mesh, d[j], ax=ax[0,j], logScale=True, cMap='jet')
            ax[0,j].set_title(tit + f'It.{j*step+1}')
            for i in range(n):    
                pg.show(self.mesh, self.weightsHistory[i][j], ax=ax[i+1,j], cMap='jet', logScale=True)
                ax[i+1,j].set_title(f'{self.names[i]} weights It.{j*step+1}')
         