module Pade
using LinearAlgebra
export PadeFermi

function B(m)
        return 2.0*m-1
end

function PadeFermi(npade)
        diagon = fill(0.0,2*npade)
        ofdiag = [1/( sqrt( B(i)*B(i+1) ) ) for i in 1:2*npade-1]
        p = SymTridiagonal(diagon,ofdiag)
        ewp,evp = eigen(p)
        diagon = ewp
        eig_denum = diagon[npade+1:2*npade]
        eig_denum = 2.0 ./eig_denum
        ofdiag = [1/( sqrt( B(i)*B(i+1) ) ) for i in 2:2*npade-1]
        diagon = fill(0.0,2*npade-1)
        p = SymTridiagonal(diagon,ofdiag)
        ewp,evp = eigen(p)
        diagon = ewp
        ofdiag = fill(0.0,2*npade-2)
        eig_num = diagon[npade+1:2*npade-1]
        eig_num = 2.0 ./ eig_num
        kappa = [1.0 for i in 1:npade]
        for l1 in 1:npade
                for l2 in 1:npade-1
                        if l2 != l1
                                kappa[l1] = kappa[l1] * ( ( eig_num[l2]^2 - eig_denum[l1]^2 ) / ( eig_denum[l2]^2 - eig_denum[l1]^2 ) )
                        else
                                kappa[l1] = kappa[l1] * (eig_num[l2]^2 - eig_denum[l1]^2)
                        end
                end
                if l1 != npade
                        kappa[l1] = kappa[l1] / (eig_denum[npade]^2 - eig_denum[l1]^2)
                end
                kappa[l1] = kappa[l1] * 0.5 * npade * B(npade+1)
        end
        return kappa,eig_denum
end


end
