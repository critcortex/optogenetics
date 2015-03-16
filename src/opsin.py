from neuron import h

class Opsin:

        
    def express_opsin(self,cell,params):
        expressed_op = {}
        print '------ About to start expressing opsins ----------'
        print 'We observe: ',params['opsindict']
        for opsin in params['opsindict'].keys():
            oplist = h.List()
            print 'Expressing opsin %s --------------------'%opsin
            for subdomain in params['opsindict'][opsin].keys():
                print '- expressing %s in %s'%(opsin,subdomain)
                
                # get the idx offset. This will only be used if setting unequal illumination
                idx_offset = cell.get_idx(subdomain)[0]
                print '_idx offset = ',idx_offset
                
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
                print 'Expressed opsin %s in %g segments in subdomain %s'%(opsin,total,subdomain)
            print '... finished expressing opsin %s'%opsin
            expressed_op[opsin] = oplist
        print '------ Finished expressing opsins ------'
        return expressed_op
    
        
    
    

    def _express_opsin_section(self,segment,section,opsintype,opsindict,idx=None,cell=None,):
        """ Expresses opsin in a specific segment """
            
        opsin = getattr(h,opsintype)()
        opsin.loc(segment.x,sec=section) 
        opsin.gbar = opsindict['exp']*h.area(segment.x)
        
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
            print mid_depth, (top_depth-mid_depth), effective_illumination
            #opsin.irr_factor = effective_illumination
            
        except:
            #print('!! Required information missing for calculating effective_illumination')
            #print("!! Required fields for opdict = projection, irr_gradient, irr_surface")
            pass
            
        
        for (k,v) in opsindict.iteritems():
            # Note that we have a couple of special cases we have to ignore (and cannot delete from the dictionary)
            if k in ['exp','irr_surface','irr_gradient','projection']:
                continue
            try:
                setattr(opsin,k,v)
            except KeyError, LookupError:
                print 'Could not update for %s in opsin %s'%(k,opsintype)
                
        return opsin
        
"""
# written for when sectionlist was a SectionList    
        
    def express_opsin(self,cell,params):
        expressed_op = {}
        print '=============== About to start expressing opsins =================='
        print 'We observe: ',params['opsindict']
        for opsin in params['opsindict'].keys():
            oplist = h.List()
            print 'Expressing opsin %s --------------------'%opsin
            for subdomain in params['opsindict'][opsin].keys():
                print '- expressing %s in %s'%(opsin,subdomain)
                total = 0
                #get corresponding section list
                sectionlist = cell.domainlists[subdomain]
                for sec in sectionlist:
                    #print sec.name()
                    for (j,seg) in enumerate(sec):
                        op = self._express_opsin_section(seg,sec,opsin,params['opsindict'][opsin][subdomain])
                        oplist.append(op)
                        total +=1
                print 'Expressed opsin %s in %g segments in subdomain %s'%(opsin,total,subdomain)
            print '... finished expressing opsin %s'%opsin
            expressed_op[opsin] = oplist
        print '=============== Finished expressing opsins =================='
        return expressed_op
    
    

    def _express_opsin_section(self,segment,section,opsintype,opsindict):
        "" Expresses opsin in a specific segment ""

        opsin = getattr(h,opsintype)()
        opsin.loc(segment.x,sec=section)
        opsin.gbar = opsindict['exp']*h.area(segment.x)
                
        for (k,v) in opsindict.iteritems():
            if k=='exp':
                continue
            try:
                setattr(opsin,k,v)
            except KeyError, LookupError:
                print 'Could not update for %s in opsin %s'%(k,opsintype)
                
        return opsin
        
"""    
