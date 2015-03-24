from neuron import h
import logging

class Opsin:
    """
    
    Designed to work with corresponding .mod files within NEURON
    
    Values can be altered by supplying a dictionary of parameters containing an opsindict:
    opsindict =  {  'exp':5e-4, 
                    'irradiance' :0.005, 
                    'pulsewidth': light_dur,
                    'lightdelay':light_on,
                    'interpulse_interval':ipi,  
                    'n_pulses':n_pulses}
    
    
    """

        
    def express_opsin(self,cell,params):
        """
        Public method for expressing the opsin throughout the entire cell
        
        """
        expressed_op = {}
        logging.info('------ About to start expressing opsins ----------')
        logging.info('We observe: ',params['opsindict'])
        for opsin in params['opsindict'].keys():
            oplist = h.List()
            logging.info('Expressing opsin %s --------------------'%opsin)
            for subdomain in params['opsindict'][opsin].keys():
                logging.info('- expressing %s in %s'%(opsin,subdomain))
                
                # get the idx offset. This will only be used if setting unequal illumination
                idx_offset = cell.get_idx(subdomain)[0]
                logging.debug('_idx offset = %g'%idx_offset)
                
                total = 0  # total sections expressed so far
                #get corresponding section list
                sectionlist = cell.domainlists[subdomain]
                slcount= sectionlist.count()
                for i in range(int(slcount)):
                    #print sec.name()
                    sec = sectionlist.o(i)
                    for (j,seg) in enumerate(sec):
                        #print type(seg), type(sec)
                        op = self._express_opsin_section(seg,sec,opsin,params['opsindict'][opsin][subdomain],idx_offset+total,cell)
                        oplist.append(op)
                        total +=1
                logging.info('Expressed opsin %s in %g segments in subdomain %s'%(opsin,total,subdomain))
            logging.info('... finished expressing opsin %s'%opsin)
            expressed_op[opsin] = oplist
        logging.info('------ Finished expressing opsins ------')
        return expressed_op
    
        
    
    

    def _express_opsin_section(self,segment,section,opsintype,opsindict,idx=None,cell=None,):
        """ 
        Expresses opsin in a specific segment.  
        
        First, goes through and sets the opsin parameters. 
        Then, if there is partial illumination, goes through and calculates the effective illumination 
        seen at this segment. This is calculated in this order to avoid overwriting the calculated "effective" illumination 
        with the "true" illumination/irradiance 
        
        @param
            segment
            section
            opsintype
            opsindict
            idx
            cell
        
        """
    
        opsin = getattr(h,opsintype)()
        opsin.loc(segment.x,sec=section) 
        opsin.gbar = opsindict['exp']*h.area(segment.x)
        
        for (k,v) in opsindict.iteritems():
            # Note that we have a couple of special cases we have to ignore (and cannot delete from the dictionary)
            if k in ['exp','irr_surface','irr_gradient','projection']:
                continue
            try:
                setattr(opsin,k,v)
            except KeyError, LookupError:
                logging.warning('Could not update for %s in opsin %s'%(k,opsintype))
                
              
        try: # see if we need to effective_illumination at each depth
            projection = opsindict['projection']
            if projection is not None and cell is not None and idx is not None:
                if projection is 'z':
                    mid_depth = cell.zmid[idx] 
                    top_depth = cell.zmid.max()
                elif projection is 'y':
                    mid_depth = cell.ymid[idx] 
                    top_depth = cell.ymid.max()
                elif projection is 'x':
                    mid_depth = cell.xmid[idx]
                    top_depth = cell.xmid.max()
            
            effective_illumination = opsindict['irr_surface'] - opsindict['irr_gradient']*(top_depth-mid_depth)
            if effective_illumination <0:
                effective_illumination = 0
            logging.debug("depth debugging: %g %g %g"%( mid_depth, (top_depth-mid_depth), effective_illumination))
            opsin.irradiance = effective_illumination*opsindict['irradiance']
            
        except:
            logging.info('!! Required information missing for calculating effective_illumination')
            logging.info("!! Required fields for opdict = projection, irr_gradient, irr_surface")
  
        return opsin
        
