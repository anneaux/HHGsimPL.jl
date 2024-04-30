# ### everything to decide whether or not a given saddle point contributes.
# ### this could be implemented in various methods again. Also maybe it should give a warning if there're multiple saddle points nearby and if a Gaussian approximation is a bad idea?



# ### necklace code
# # gradient()
# mutable struct Simplex{T}
# 	const t1::T
# 	const t2::T where T <:Number
# 	active::Bool
# end



# # using LinearAlgebra

# # function Necklace(qx, Rx, Φx; ts, NN=10, ϵ=0.01, δ=0.05, Δ=0.5)
# #     gradN(q, R, Φ, a) = normalize(conj.(complex(D[I * action(q, R, Φ, a), [a1, a2]])))
    
# #     counter = 0
# #     while counter < 500
# #         counter += 1
# #         tmp = deepcopy(necklace)
# #         for i in 1:length(necklace)
# #             if real(I * action(qx, Rx, Φx, necklace[i])) < 0
# #                 necklace[i] += δ * gradN(qx, Rx, Φx, necklace[i])
# #             end
# #         end
# #         for i in 1:length(necklace)-1
# #             if norm(necklace[i] - necklace[i+1]) > Δ
# #                 insert!(necklace, i+1, (necklace[i] + necklace[i+1]) / 2)
# #             end
# #         end
# #         if necklace == tmp
# #             break
# #         end
# #     end
# #     return necklace
# # end





# ### initialise
# # function initialise_necklace(beam, NN::Int64=50)
# #     hess(q, R, Φ, a) = begin
# #         ∇²action = Complex(diff(re(action(q, R, Φ, [u1 + v1*im, u2 + v2*im]))), [u1, v1, u2, v2])
# #         hess_matrix = complex(diff(∇²action, [u1, v1, u2, v2], [u1, v1, u2, v2]))
# #         return hess_matrix
# #     end
    
# #     eigensystem = eig(hess(qx, Rx, Φx, ts))
# #     sorted_eigensystem = eigensystem.vectors[:, sortperm(real.(eigensystem.values))]
# #     for i in 1:length(sorted_eigensystem)
# #         sorted_eigensystem[:, i] = toComplex(normalize(sorted_eigensystem[:, i]))
# #     end
    
# #     necklace = [ts .+ ϵ * (cos(θ) * sorted_eigensystem[:, 3] + sin(θ) * sorted_eigensystem[:, 4]) for θ in range(0, stop=2π, length=NN+1)]

# #     return necklace
# # end



# ### subdevide



# ### flow



# ### clean