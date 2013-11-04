from neuron import h

class Opsin:

        
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
        """ Expresses opsin in a specific segment """

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
        
    
